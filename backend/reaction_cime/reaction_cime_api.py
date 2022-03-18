import pickle
from flask import Blueprint, request, current_app, abort, jsonify, Response, stream_with_context
import logging
from flask.helpers import make_response, send_file
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from .helper_functions import (circ_radius_to_radius, preprocess_dataset, get_mcs, smiles_to_base64, aggregate_by_col_interpolate, 
                                aggregate_by_col, rescale_and_encode, hex_aggregate_by_col, isInsideHex, circ_radius_to_radius)
import json

_log = logging.getLogger(__name__)

reaction_cime_api = Blueprint('reaction_cime', __name__)

@reaction_cime_api.route('/hello', methods=['GET'])
def hello():
    return "Hello World"


# --------- server management --------- 
from flask import render_template
@reaction_cime_api.route("/")
def index():
    return render_template("jku-vds-lab/reaction-cime/index.html")

@reaction_cime_api.route('/index.html')
def index_2():
    return render_template('jku-vds-lab/reaction-cime/index.html')


# @reaction_cime_api.route('/favicon.ico')
# def favicon():
#     return render_template('jku-vds-lab/reaction-cime/favicon.ico')

# --------------- database ------------------

@reaction_cime_api.route('/get_uploaded_files_list', methods=['GET'])
def get_uploaded_files_list():
    return jsonify(get_cime_dbo().get_table_names())

@reaction_cime_api.route('/delete_file/<filename>', methods=['GET'])
def delete_uploaded_file(filename):
    deleted = get_cime_dbo().drop_table(filename)
    print("----")
    print(deleted)
    # TODO: Refactor to true and false instead of strings
    return {"deleted": "true" if deleted else "false"}


def get_cime_dbo():
    return current_app.config['REACTION_CIME_DBO']


# --------------- process dataset ------------------
import os
import time
import numpy as np
@reaction_cime_api.route('/upload_csv', methods=['OPTIONS', 'POST'])
def upload_csv():
    start_time = time.time()
    _log.info(f'Received new csv to upload')
    fileUpload = request.files.get("myFile")

    # --- check if file is corrupt
    if not fileUpload or fileUpload.filename == '':
        abort('No valid file provided')
    fileUpload.seek(0, 2) # sets file's current position at 0th offset from the end (2)
    file_length = fileUpload.tell() # get absolute offset of current position
    supposed_file_size = int(request.form.get("file_size"))
    if file_length != supposed_file_size: # the uploaded file does not correspond to the original file
        # abort("there was a problem with the file upload. please try again") # not working?
        return {"error": "there was a problem with the file upload. please try again"}
    fileUpload.seek(0, 0) # resets file's current position at 0th offset from start (0)
    # ---

    _log.info(f'Uploading file {fileUpload.filename}')
    df = pd.read_csv(request.files.get('myFile'))

    if "x" not in df.columns or "y" not in df.columns:
        print("--- randomly init x and y coordinates")
        df["x"] = np.random.uniform(-50, 50, len(df))
        df["y"] = np.random.uniform(-50, 50, len(df))


    print("--- save file to database")
    filename = fileUpload.filename
    filename = "_".join(filename.split(".")[0:-1])
    if filename in get_cime_dbo().get_table_names():
        filename = "%s%i"%(filename, time.time())

    get_cime_dbo().save_dataframe(df, filename)
    
    delta_time = time.time()-start_time
    print("--- took", time.strftime('%H:%M:%S', time.gmtime(delta_time)), "to upload file %s"%filename)
    print("--- took %i min %f s to upload file %s"%(delta_time/60, delta_time%60, filename))

    return {
        "filename": fileUpload.filename,
        "id": filename
    }



import pandas as pd
from io import StringIO
@reaction_cime_api.route('/get_poi_csv/<filename>', methods=['GET'])
def get_points_of_interest(filename):
    # domain = pd.read_csv("./temp-files/%s"%filename)
    # domain = get_cime_dbo().get_dataframe_from_table(filename)

    poi_domain = get_poi_df_from_db(filename, get_cime_dbo())

    poi_domain = preprocess_dataset(poi_domain)

    csv_buffer = StringIO()
    poi_domain.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()



def get_poi_df_from_db(filename, cime_dbo):
    # TODO: make dynamic, by which feature we want to filter (e.g. user could change the settings in the front-end maybe with parallel coordinates?)
    # example code for distance filter
    poi_domain = cime_dbo.get_dataframe_from_table_filter(filename, "yield > 0")
    return poi_domain

def get_poi_mask(filename, cime_dbo):
    mask = cime_dbo.get_filter_mask(filename, "yield > 0")
    return mask

@reaction_cime_api.route('/get_k_nearest_from_csv/<filename>/<x>/<y>/<k>', methods=['GET'])
def get_k_nearest_points(filename, x, y, k):
    """
    k nearest points to the point x,y in the database filename will be returned using a StringIO buffer.
    This is implemented via an SQLite query to the database.

    Parameters
    ----------
    filename:
        link to the dataset.

    x:
        The x coordinate of the point for which the nearest neighbours should be returned
    y:
        The y coordinate of the point for which the nearest neighbours should be returned
    k:
        The amount of neighbours to return

    Returns
    ----------
    str
        csv created from the filtered database entries
    """
    
    # calculates squared euclidean distance, orders the db table by this distance, and returns k first entries
    nearest_points_domain = get_cime_dbo().get_dataframe_from_table_complete_filter(filename, 'ORDER BY ((x-('+str(x)+'))*(x-('+str(x)+')))+((y-('+str(y)+'))*(y-('+str(y)+'))) LIMIT '+str(k)+';')

    response = make_response(nearest_points_domain.to_csv(index=False))
    response.headers["Content-Disposition"] = "attachment; filename=%s_kNN_%s_%s.csv"%(filename, x, y)
    response.headers["Content-type"] = "text/csv"
    return response

@reaction_cime_api.route('/get_radius_from_csv/<filename>/<x>/<y>/<d>', methods=['GET'])
def get_points_given_radius(filename, x, y, r):
    """
    Returns all points in the database filename that have distance no greater than r to the point at coordinates (x,y).

    Parameters
    ----------
    filename:
        link to the dataset.

    x:
        x coordinate of the point for which the nearest neighbours should be returned
    y:
        y coordinate of the point for which the nearest neighbours should be returned
    r:
        radius around the specified coordinate in which to return data points

    Returns
    ----------
    str
        csv created from the filtered database entries
    """
    # squared radius for SQLite query filter condition
    r2 = int(r) * int(r)
    # filters using euclidean distance (squared on both sides since SQLite does not have extended maths enabled for SQRT)
    nearest_points_domain = get_cime_dbo().get_dataframe_from_table_filter(filename, '((x-('+str(x)+'))*(x-('+str(x)+')))+((y-('+str(y)+'))*(y-('+str(y)+'))) < '+str(r2))

    response = make_response(nearest_points_domain.to_csv(index=False))
    response.headers["Content-Disposition"] = "attachment; filename=%s_kNN_%s_%s.csv"%(filename, x, y)
    response.headers["Content-type"] = "text/csv"
    return response


# deprecated
@reaction_cime_api.route('/get_agg_csv/<filename>/<col_name>', methods=['GET'])
def get_aggregated_dataset(filename, col_name):

    range = {"x_min": request.args.get("x_min"), 
        "x_max": request.args.get("x_max"),
        "y_min": request.args.get("y_min"), 
        "y_max": request.args.get("y_max"),
        }
    filter = "x > {x_min} and x < {x_max} and y > {y_min} and y < {y_max}".format(**range)
    agg_domain = get_cime_dbo().get_dataframe_from_table_filter(filename, filter, columns=["x", "y", col_name])
    agg_df = aggregate_by_col_interpolate(agg_domain, col_name, sample_size=200) # TODO: dynamic sample_size

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


dataset_cache = {} # structure: {"filename": dataset} # if "retrieve_cols" are contained in dataset.columns then the dataset is returned, otherwise dataset with all "cache_cols" are loaded into the cache
def handle_dataset_cache(filename, cols=[]):
    if filename not in dataset_cache.keys():
        dataset_cache[filename] = None

    all_cols = list(set(cols + ["x", "y"]))
    # dataset is not cached and has to be loaded
    if dataset_cache[filename] is None or not set(all_cols).issubset(set(dataset_cache[filename].columns)):
        print("---dataset is not cached and has to be loaded to memory")
        print(all_cols)
        if(dataset_cache[filename] is not None): # add currently cached columns again to cache
            all_cols = list(set(all_cols + list(dataset_cache[filename].columns)))
        print(all_cols)
        dataset_cache[filename] = get_cime_dbo().get_dataframe_from_table(filename, columns=all_cols)
    else:
        print("---dataset cached!")
    # return cached version of dataset
    return dataset_cache[filename]

def reset_dataset_cache(filename=None):
    global dataset_cache
    if filename is None:
        dataset_cache = {}
    else:
        dataset_cache[filename] = None


@reaction_cime_api.route('/update_cache/<filename>', methods=['GET'])
def update_dataset_cache(filename):
    cache_cols = request.args.getlist("cache_cols")
    handle_dataset_cache(filename, cache_cols)
    return {"msg": "ok"}

@reaction_cime_api.route('/get_agg_csv_cached/<filename>', methods=['GET'])
def get_aggregated_dataset_cached(filename):
    # value_col_name = request.args.get("value_col")
    # uncertainty_col_name = request.args.get("uncertainty_col")
    retrieve_cols = request.args.getlist("retrieve_cols")
    cache_cols = request.args.getlist("cache_cols") # if value col or uncertainty col are not cached, we use this list of columns to prepare the new cached dataset
    sample_size = request.args.get("sample_size", default=200, type=int)
    range = {"x_min": float(request.args.get("x_min")),
        "x_max": float(request.args.get("x_max")),
        "y_min": float(request.args.get("y_min")),
        "y_max": float(request.args.get("y_max")),
        }

    agg_domain = handle_dataset_cache(filename, retrieve_cols + cache_cols)
    agg_domain = agg_domain[(agg_domain["x"] < range["x_max"]) * (agg_domain["x"] > range["x_min"]) * (agg_domain["y"] < range["y_max"]) * (agg_domain["y"] > range["y_min"])]

    agg_df = aggregate_by_col_interpolate(agg_domain, retrieve_cols, sample_size=sample_size)
    # agg_df = aggregate_by_col(agg_domain, retrieve_cols, sample_size=sample_size)

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


@reaction_cime_api.route('/get_hex_agg/<filename>', methods=['GET'])
def get_hexagonal_aggregation(filename):
    retrieve_cols = request.args.getlist("retrieve_cols")
    aggregation_methods = request.args.getlist("aggregation_methods")
    cache_cols = request.args.getlist("cache_cols") # if value col or uncertainty col are not cached, we use this list of columns to prepare the new cached dataset
    sample_size = 20#request.args.get("sample_size", default=20, type=int)
    range = {"x_min": float(request.args.get("x_min")),
        "x_max": float(request.args.get("x_max")),
        "y_min": float(request.args.get("y_min")),
        "y_max": float(request.args.get("y_max")),
        }
    
    agg_domain = handle_dataset_cache(filename, retrieve_cols + cache_cols)
    agg_domain = agg_domain[(agg_domain["x"] < range["x_max"]) * (agg_domain["x"] > range["x_min"]) * (agg_domain["y"] < range["y_max"]) * (agg_domain["y"] > range["y_min"])]

    agg_df, wrong_points = hex_aggregate_by_col(agg_domain, retrieve_cols, aggregation_methods, range=None, sample_size=sample_size) # TODO: what works better: range set to the boundaries of the dataset, or range set to the boundaries of the screen i.e. range=range

    wrong_df = wrong_points[list(set(retrieve_cols + ["x", "y"]))]
    wrong_df["hex"] = False

    agg_df["hex"] = True
    
    agg_df = agg_df.append(wrong_df, ignore_index=True)

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


@reaction_cime_api.route('/get_value_range/<filename>/<col_name>', methods=['GET'])
def get_value_range(filename, col_name):
    range = get_cime_dbo().get_value_range_from_table(filename, col_name).iloc[0]
    return {"min": float(range["min"]), "max": float(range["max"])}


@reaction_cime_api.route('/get_category_count/<filename>/<col_name>', methods=["GET"])
def get_category_count(filename, col_name):
    # cat_count = get_cime_dbo().get_category_count(filename, col_name)
    # return json.dumps(cat_count.to_dict('records')).encode('utf-8')
    df = handle_dataset_cache(filename, [col_name])
    df = df[[col_name]]
    df["count"] = 0 # dummy column, such that it can be aggregated by count
    cat_count = df.groupby(col_name, as_index=False).count()
    return json.dumps(cat_count.to_dict('records')).encode('utf-8')

@reaction_cime_api.route('/get_category_count_of_hex/<filename>/<col_name>', methods=["GET"])
def get_category_count_of_hex(filename, col_name):
    hex_x = float(request.args.get("x"))
    hex_y = float(request.args.get("y"))
    hex_circ_radius = float(request.args.get("circ_radius"))

    # df = get_cime_dbo().get_dataframe_from_table(filename, columns=["x", "y", col_name])
    df = handle_dataset_cache(filename, ["x", "y", col_name])
    
    window = isInsideHex(df["x"], df["y"], hex_x=hex_x, hex_y=hex_y, radius=circ_radius_to_radius(hex_circ_radius), circ_radius=hex_circ_radius)
    df = df[window][[col_name]]
    df["count"] = 0 # dummy column, such that it can be aggregated by count
    cat_count = df.groupby(col_name, as_index=False).count()
    return json.dumps(cat_count.to_dict('records')).encode('utf-8')

@reaction_cime_api.route('/get_density/<filename>/<col_name>', methods=["GET"])
def get_density(filename, col_name):
    # data = get_cime_dbo().get_dataframe_from_table(filename, columns=[col_name])[col_name]
    data = handle_dataset_cache(filename, [col_name])[col_name]

    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    x_vals = np.linspace(min(data), max(data), 100)
    y_vals = density(x_vals)

    return {"x_vals": list(x_vals), "y_vals": list(y_vals)}

@reaction_cime_api.route('/get_density_of_hex/<filename>/<col_name>', methods=["GET"])
def get_density_of_hex(filename, col_name):
    hex_x = float(request.args.get("x"))
    hex_y = float(request.args.get("y"))
    hex_circ_radius = float(request.args.get("circ_radius"))

    # df = get_cime_dbo().get_dataframe_from_table(filename, columns=["x", "y", col_name])
    df = handle_dataset_cache(filename, ["x", "y", col_name])
    
    window = isInsideHex(df["x"], df["y"], hex_x=hex_x, hex_y=hex_y, radius=circ_radius_to_radius(hex_circ_radius), circ_radius=hex_circ_radius)
    data = df[window][col_name]

    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    x_vals = np.linspace(min(data), max(data), 100)
    y_vals = density(x_vals)

    return {"x_vals": list(x_vals), "y_vals": list(y_vals)}


from openTSNE import TSNE
import umap
from sklearn.decomposition import PCA
import gower
# depricated
@reaction_cime_api.route('/project_dataset', methods=['OPTIONS', 'POST'])
def project_dataset():
    filename = request.form.get("filename")
    params = json.loads(request.form.get("params"))
    selected_feature_info = json.loads(request.form.get("selected_feature_info"))
    proj_df = get_cime_dbo().get_dataframe_from_table(filename, columns=list(selected_feature_info.keys()))

    # TODO: if params["useSelection"]: only project selected points --> do this in front-end? do this at all?

    # handle custom initialization of coordinates
    initialization = None
    if params["seeded"]: # take custom seed from current coordinates
        initialization = get_poi_df_from_db(filename, get_cime_dbo())[["x","y"]].values


    # rescale numerical values and encode categorical values
    proj_df, categorical_features = rescale_and_encode(proj_df, params, selected_feature_info)


    # handle custom metrics
    metric = params["distanceMetric"]
    if metric == None or metric == "":
        metric = "euclidean"
    if metric == "gower": # we precompute the similarity matrix with the gower metric
        metric = "precomputed"
        normalized_values = gower.gower_matrix(proj_df, cat_features=categorical_features)
    else: # otherwise, the similarity can be done by a pre-configured function and we just hand over the values
        normalized_values = proj_df.values



    # project the data
    if params["embedding_method"] == "umap":
        if initialization == None:
            initialization = "spectral"

        proj = umap.UMAP(n_neighbors=int(params["nNeighbors"]), n_components=2, metric=metric, n_epochs=int(params["iterations"]), init=initialization, verbose=True) # output_metric="euclidean" --> ? learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option 
        proj_data = proj.fit_transform(normalized_values)
    elif params["embedding_method"] == "tsne":
        if initialization == None:
            initialization = "pca"
        proj = TSNE(2, metric=metric, perplexity=int(params["perplexity"]), n_iter=int(params["iterations"]), initialization=initialization, verbose=True) #, learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option 
        proj_data = proj.fit(normalized_values)
    else:
        print("--- defaulting to PCA embedding")
        pca = PCA(n_components=2)
        proj_data = pca.fit_transform(normalized_values)


    # update the coordinates in the dataset
    start_time = time.time()
    get_cime_dbo().update_row_bulk(filename, proj_df.index, {"x":proj_data[:,0], "y": proj_data[:,1]})
    delta_time = time.time()-start_time
    print("--- took", time.strftime('%H:%M:%S', time.gmtime(delta_time)), "to update database")
    print("--- took %i min %f s to update database"%(delta_time/60, delta_time%60))

    print("return poi positions")
    return {"done": "true", "embedding": get_poi_df_from_db(filename, get_cime_dbo())[["x","y"]].to_dict('records')}



@reaction_cime_api.route('/terminate_projection_thread/<filename>', methods=["GET"])
def terminate_projection_thread(filename):
    current_app.config["TERMINATE_PROJECTION_" + filename] = True
    return {"msg": "ok"}

    
import threading
import ctypes
import time
  
# https://www.geeksforgeeks.org/python-different-ways-to-kill-a-thread/
class ProjectionThread(threading.Thread):
    def __init__(self, filename, params, selected_feature_info, cime_dbo):
        threading.Thread.__init__(self)
        self.filename = filename
        self.params = params
        self.selected_feature_info = selected_feature_info
        self.cime_dbo = cime_dbo

        self.current_step = "0"
        self.done = False
        self.emb = None
        self.msg = "init..."
             
    def run(self):
        try:

            if self.params["embedding_method"] == "rmOverlap":
                self.msg = "load dataset..."
                proj_df = self.cime_dbo.get_dataframe_from_table(self.filename, columns=["x", "y"])

                self.msg = "calc embedding..."
                from .dgrid import DGrid
                # TODO: make parameters dynamic
                icon_width = 1
                icon_height = 1
                delta = 6
                start_time = time.time()
                proj_data = DGrid(icon_width=icon_width, icon_height=icon_height, delta=delta).fit_transform(proj_df[["x","y"]].values)
                delta_time = time.time()-start_time
                print("--- took", time.strftime('%H:%M:%S', time.gmtime(delta_time)), "calculate remove overlap")
            else:
                self.msg = "load dataset..."
                proj_df = self.cime_dbo.get_dataframe_from_table(self.filename, columns=list(self.selected_feature_info.keys()))

                self.msg = "seeding..."
                # handle custom initialization of coordinates
                initialization = None
                if self.params["seeded"]: # take custom seed from current coordinates
                    # initialization = get_poi_df_from_db(self.filename, self.cime_dbo)[["x","y"]].values # does this make sense?
                    initialization = proj_df[["x","y"]].values

                # TODO: if params["useSelection"]: only project selected points --> do this in front-end? do this at all?

                self.msg = "rescale and encode..."
                # rescale numerical values and encode categorical values
                proj_df, categorical_features = rescale_and_encode(proj_df, self.params, self.selected_feature_info)


                self.msg = "calc metric..."
                # handle custom metrics
                metric = self.params["distanceMetric"]
                if metric == None or metric == "":
                    metric = "euclidean"

                if metric == "gower": # we precompute the similarity matrix with the gower metric
                    metric = "precomputed"
                    normalized_values = gower.gower_matrix(proj_df, cat_features=categorical_features)

                    if initialization == "pca" and self.params["embedding_method"] == "tsne":
                        print("--- calculating PCA for tSNE with Gowers distance")
                        pca = PCA(n_components=2)
                        initialization = pca.fit_transform(proj_df.values)
                    
                else: # otherwise, the similarity can be done by a pre-configured function and we just hand over the values
                    normalized_values = proj_df.values

                
                self.msg = "calc embedding..."
                # project the data
                # --- umap
                if self.params["embedding_method"] == "umap":
                    if initialization == None:
                        initialization = "spectral"

                    # https://pypi.org/project/tqdm/#parameters
                    # use TQDM to update steps for front-end
                    class CustomTQDM:
                        def __init__(self, parent):
                            self.parent = parent

                        def write(self, str):
                            try:
                                int(str) # if it is not possible to parse it as int, we reached the end?
                                self.parent.current_step = str
                            except ValueError:
                                print(str)
                            
                        def flush(self):
                            pass

                    proj = umap.UMAP(n_neighbors=int(self.params["nNeighbors"]), n_components=2, metric=metric, n_epochs=int(self.params["iterations"]), init=initialization, tqdm_kwds={"bar_format": '{n_fmt}' ,"file": CustomTQDM(self)}, verbose=True) # output_metric="euclidean" --> ? learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option 
                    proj_data = proj.fit_transform(normalized_values)

                # --- tsne
                elif self.params["embedding_method"] == "tsne":
                    # https://github.com/pavlin-policar/openTSNE/blob/master/openTSNE/callbacks.py
                    from openTSNE.callbacks import Callback
                    class CustomCallback(Callback):
                        def __init__(self, parent):
                            self.parent = parent
                            
                        def optimization_about_to_start(self):
                            """This is called at the beginning of the optimization procedure."""
                            print("start optimization")
                            pass
                            
                        def __call__(self, iteration, error, embedding): # error is current KL divergence of the given embedding
                            poi_mask = get_poi_mask(self.parent.filename, self.parent.cime_dbo)
                            self.parent.current_step = iteration
                            self.parent.emb = [{"x": row[0], "y": row[1]} for row in embedding[poi_mask["mask"]]]
                    #         return True # with this optimization will be interrupted

                    if initialization == None:
                        initialization = "pca"
                    proj = TSNE(2, metric=metric, callbacks_every_iters=10, callbacks=[CustomCallback(self)], perplexity=int(self.params["perplexity"]), n_iter=int(self.params["iterations"]), initialization=initialization, verbose=True) #, learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option 
                    proj_data = proj.fit(normalized_values)

                # --- pca
                else:
                    print("--- defaulting to PCA embedding")
                    pca = PCA(n_components=2)
                    proj_data = pca.fit_transform(normalized_values)
            
            self.current_step = self.params["iterations"]
            self.msg = "save embedding..."
            # update the coordinates in the dataset
            start_time = time.time()
            self.cime_dbo.update_row_bulk(self.filename, proj_df.index, {"x":proj_data[:,0], "y": proj_data[:,1]})
            reset_dataset_cache(self.filename) # reset cached data
            delta_time = time.time()-start_time
            print("--- took", time.strftime('%H:%M:%S', time.gmtime(delta_time)), "to update database")
            print("--- took %i min %f s to update database"%(delta_time/60, delta_time%60))

            self.msg = "finished!"

        except Exception as e:
            self.msg = "error during projection"
            print(e)

        finally:
            self.done = True
    
    def get_id(self):
 
        # returns id of the respective thread
        if hasattr(self, '_thread_id'):
            return self._thread_id
        for id, thread in threading._active.items():
            if thread is self:
                return id
  
    def raise_exception(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id,
              ctypes.py_object(SystemExit))
        if res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
            print('Exception raise failure')

        self.done = True
        self.msg = "client-side termination"



# differnt approach could have been with sockets: https://www.shanelynn.ie/asynchronous-updates-to-a-webpage-with-flask-and-socket-io/
@reaction_cime_api.route('/project_dataset_async', methods=['OPTIONS', 'POST'])
def project_dataset_async():
    filename = request.form.get("filename")
    params = json.loads(request.form.get("params"))
    selected_feature_info = json.loads(request.form.get("selected_feature_info"))
    cime_dbo = get_cime_dbo()

    current_app.config["TERMINATE_PROJECTION_" + filename] = False

    def my_generator():
        
        proj = ProjectionThread(filename, params, selected_feature_info, cime_dbo)
        proj.start()
        while True:
            time.sleep(2)
            terminate = False
            if "TERMINATE_PROJECTION_" + filename in current_app.config.keys():
                terminate = current_app.config["TERMINATE_PROJECTION_" + filename]

            # https://stackoverflow.com/questions/41146144/how-to-fix-assertionerror-value-must-be-bytes-error-in-python2-7-with-apache
            yield json.dumps({"step": proj.current_step, "msg": proj.msg, "emb": proj.emb}).encode('utf-8')
            if terminate:
                proj.raise_exception()
                proj.join()
            if proj.done: #proj.current_step >= params["iterations"] or proj.interrupt:
                yield json.dumps({"step": None, "msg": None, "emb": get_poi_df_from_db(filename, cime_dbo)[["x","y"]].to_dict('records')}).encode('utf-8')
                break
    
    # https://flask.palletsprojects.com/en/1.1.x/patterns/streaming/
    return Response(stream_with_context(my_generator()))
    




# --------------- chem specific ------------------- 

@reaction_cime_api.route('/get_mol_img', methods=['OPTIONS', 'POST'])
def smiles_to_img_post():
    if request.method == 'POST':
        smiles = request.form.get("smiles")
        img = smiles_to_base64(smiles)
        return {"data": img}
    else:
        return {}

@reaction_cime_api.route('/get_common_mol_img', methods=['OPTIONS', 'POST'])
def smiles_list_to_common_substructure_img():
    if request.method == 'POST':
        smiles_list = request.form.getlist("smiles_list")
        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}
        if len(smiles_list) == 1:
            ret = smiles_to_base64(smiles_list[0])
            return {"data": ret, "smiles": smiles_list[0]}

        mol_lst = []
        error_smiles = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_lst.append(mol)
            else:
                error_smiles.append(smiles)

        m = get_mcs(mol_lst)
        pil_img = Draw.MolToImage(m)

        buffered = BytesIO()
        pil_img.save(buffered, format="JPEG")
        img_str = base64.b64encode(buffered.getvalue())
        return {"data": img_str.decode("utf-8"), "smiles": Chem.MolToSmiles(m)}
    else:
        return {}


@reaction_cime_api.route('/get_substructure_count', methods=['OPTIONS', 'POST'])
def smiles_list_to_substructure_count():
    if request.method == 'POST':
        smiles_list = request.form.get("smiles_list").split(",")
        filter_smiles = request.form.get("filter_smiles")

        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}

        patt = Chem.MolFromSmiles(filter_smiles)
        if patt:
            substructure_counts = [(smiles, len(Chem.MolFromSmiles(smiles).GetSubstructMatch(
                patt))) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
            return {"substructure_counts": substructure_counts}
        return {"error": "invalid SMILES filter"}
    else:
        return {}

