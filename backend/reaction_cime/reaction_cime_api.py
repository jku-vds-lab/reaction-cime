import pickle
from flask import Blueprint, request, current_app, abort, jsonify, Response, stream_with_context
import logging
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from io import BytesIO
import base64
from .helper_functions import generate_rename_list, get_mcs, smiles_to_base64, aggregate_col, rescale_and_encode
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

@reaction_cime_api.route('/manifest.json')
def manifest():
    return render_template('jku-vds-lab/reaction-cime/manifest.json')

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

    print(filename)
    poi_domain = get_poi_df_from_db(filename, get_cime_dbo())
    print(len(poi_domain))

    new_cols = generate_rename_list(poi_domain)
    poi_domain.columns = new_cols

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
        StringIO.getValue() of the buffered csv created from the filtered database entries
    """
    
    # calculates squared euclidean distance, orders the db table by this distance, and returns k first entries
    nearest_points_domain = get_cime_dbo().get_dataframe_from_table_filter(filename, 'ORDER BY ((x-('+str(x)+'))*(x-('+str(x)+')))+((y-('+str(y)+'))*(y-('+str(y)+'))) LIMIT '+str(k)+';', where=False)

    csv_buffer = StringIO()
    nearest_points_domain.to_csv(csv_buffer, index=False)
    
    return csv_buffer.getvalue()

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
        StringIO.getValue() of the buffered csv created from the filtered database entries
    """
    # squared radius for SQLite query filter condition
    r2 = int(r) * int(r)
    # filters using euclidean distance (squared on both sides since SQLite does not have extended maths enabled for SQRT)
    nearest_points_domain = get_cime_dbo().get_dataframe_from_table_filter(filename, '((x-('+str(x)+'))*(x-('+str(x)+')))+((y-('+str(y)+'))*(y-('+str(y)+'))) < '+str(r2))

    csv_buffer = StringIO()
    nearest_points_domain.to_csv(csv_buffer, index=False)
    
    return csv_buffer.getvalue()

@reaction_cime_api.route('/get_agg_csv/<filename>/<col_name>', methods=['GET'])
def get_aggregated_dataset(filename, col_name):
    agg_domain = get_cime_dbo().get_dataframe_from_table(filename, columns=["x", "y", col_name])
    agg_df = aggregate_col(agg_domain, col_name, sample_size=200) # TODO: dynamic sample_size

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


from openTSNE import TSNE
import umap
from sklearn.decomposition import PCA
import gower
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

            self.msg = "load dataset..."
            self.proj_df = self.cime_dbo.get_dataframe_from_table(self.filename, columns=list(self.selected_feature_info.keys()))

            self.msg = "seeding..."
            # handle custom initialization of coordinates
            initialization = None
            if self.params["seeded"]: # take custom seed from current coordinates
                initialization = get_poi_df_from_db(self.filename, self.cime_dbo)[["x","y"]].values

            # TODO: if params["useSelection"]: only project selected points --> do this in front-end? do this at all?

            self.msg = "rescale and encode..."
            # rescale numerical values and encode categorical values
            proj_df, categorical_features = rescale_and_encode(self.proj_df, self.params, self.selected_feature_info)


            self.msg = "calc metric..."
            # handle custom metrics
            metric = self.params["distanceMetric"]
            if metric == None or metric == "":
                metric = "euclidean"
            if metric == "gower": # we precompute the similarity matrix with the gower metric
                metric = "precomputed"
                normalized_values = gower.gower_matrix(proj_df, cat_features=categorical_features)
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

            self.msg = "save embedding..."
            # update the coordinates in the dataset
            start_time = time.time()
            self.cime_dbo.update_row_bulk(self.filename, proj_df.index, {"x":proj_data[:,0], "y": proj_data[:,1]})
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
    

# @reaction_cime_api.route('/project_dataset_async_test', methods=['OPTIONS', 'POST'])
# def project_dataset_async_test():

#     def my_generator():
        
#         yield json.dumps({"step":1, "emb":"hello"}).encode('utf-8')
#         time.sleep(5)
#         yield json.dumps({"step":2, "emb":"hello2"}).encode('utf-8')
#         yield json.dumps({"step":3, "emb":"hello3"}).encode('utf-8')
#         time.sleep(1)
#         yield json.dumps({"step":5, "emb":"hello4"}).encode('utf-8')
#     return Response(my_generator())


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

