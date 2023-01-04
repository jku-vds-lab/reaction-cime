import base64
import ctypes
import json
import logging
import os
import threading
import time
from io import BytesIO, StringIO
from typing import Optional

import gower
import hdbscan
import numpy as np
import pandas as pd
import umap
from flask import Blueprint, abort, current_app, jsonify, request, stream_with_context
from flask.helpers import make_response, send_file
from flask.wrappers import Response
from openTSNE import TSNE
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.decomposition import PCA

from .helper_functions import (
    aggregate_by_col_interpolate,
    circ_radius_to_radius,
    cycle_column,
    get_mcs,
    hex_aggregate_by_col,
    is_inside_hex,
    preprocess_dataset,
    rescale_and_encode,
    smiles_to_base64,
)

_log = logging.getLogger(__name__)

reaction_cime_api = Blueprint("reaction_cime", __name__)


@reaction_cime_api.route("/hello", methods=["GET"])
def hello():
    return "Hello World"


# region --------------- database ---------------


def get_cime_dbo():
    return current_app.config["REACTION_CIME_DBO"]


@reaction_cime_api.route("/get_uploaded_files_list", methods=["GET"])
def get_uploaded_files_list():
    return jsonify(get_cime_dbo().get_table_names())


@reaction_cime_api.route("/delete_file/<filename>", methods=["GET"])
def delete_uploaded_file(filename):
    deleted = get_cime_dbo().drop_table(filename)
    if deleted:
        delete_poi_exceptions(filename)
        delete_poi_constraints(filename)
    print("----")
    print(deleted)
    # TODO: Refactor to true and false instead of strings
    return {"deleted": "true" if deleted else "false"}


# endregion

# region --------------- caching ---------------
dataset_cache = (
    {}
)  # structure: {"filename": dataset} # if "retrieve_cols" are contained in dataset.columns then the dataset is returned, otherwise dataset with all "cache_cols" are loaded into the cache


def handle_dataset_cache(filename, cols=[], x_channel="x", y_channel="y"):
    if filename not in dataset_cache.keys():
        dataset_cache[filename] = None

    all_cols = list(set(cols + [x_channel, y_channel, "x", "y"]))
    # dataset is not cached and has to be loaded
    if dataset_cache[filename] is None or not set(all_cols).issubset(set(dataset_cache[filename].columns)):
        print("---dataset is not cached and has to be loaded to memory")
        print(all_cols)
        if dataset_cache[filename] is not None:  # add currently cached columns again to cache
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


@reaction_cime_api.route("/update_cache/<filename>", methods=["GET"])
def update_dataset_cache(filename):
    cache_cols = request.args.getlist("cache_cols")
    handle_dataset_cache(filename, cache_cols)
    return {"msg": "ok"}


# endregion

# region --------------- process dataset ------------------


@reaction_cime_api.route("/upload_csv", methods=["OPTIONS", "POST"])
def upload_csv():
    start_time = time.time()
    _log.info("Received new csv to upload")
    file_upload = request.files.get("myFile")

    # --- check if file is corrupt
    if not file_upload or file_upload.filename == "":
        abort(500, "No valid file provided")
    file_upload.seek(0, 2)  # sets file's current position at 0th offset from the end (2)
    file_length = file_upload.tell()  # get absolute offset of current position
    supposed_file_size = int(request.form["file_size"])
    if file_length != supposed_file_size:  # the uploaded file does not correspond to the original file
        # abort("there was a problem with the file upload. please try again") # not working?
        return {"error": "there was a problem with the file upload. please try again"}
    file_upload.seek(0, 0)  # resets file's current position at 0th offset from start (0)
    # ---

    _log.info(f"Uploading file {file_upload.filename}")
    df = pd.read_csv(request.files["myFile"].stream)

    if "x" not in df.columns or "y" not in df.columns:
        print("--- randomly init x and y coordinates")
        df["x"] = np.random.uniform(-50, 50, len(df))
        df["y"] = np.random.uniform(-50, 50, len(df))

    print("--- save file to database")
    filename = file_upload.filename or ""
    filename = "_".join(filename.split(".")[0:-1])
    if filename in get_cime_dbo().get_table_names():
        filename = "%s%i" % (filename, time.time())

    get_cime_dbo().save_dataframe(df, filename)

    # create default constraints file
    save_poi_constraints(filename)
    save_poi_exceptions(filename)

    delta_time = time.time() - start_time
    print("--- took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), "to upload file %s" % filename)
    print("--- took %i min %f s to upload file %s" % (delta_time / 60, delta_time % 60, filename))

    return {
        "filename": file_upload.filename,
        "id": filename,
    }


@reaction_cime_api.route("/get_no_datapoints/<filename>", methods=["GET"])
def get_no_datapoints(filename):
    count = get_cime_dbo().get_no_points_from_table(filename)["count"][0]
    return {"no_datapoints": int(count)}


# endregion

# region --------------- Points of Interest ------------------

MAX_POINTS = 10000


def delete_poi_exceptions(filename):
    file = f'{current_app.config["tempdir"]}{filename}_exceptions.csv'
    if os.path.exists(file):
        os.remove(file)
    else:
        print("The file does not exist")


def delete_poi_constraints(filename):
    file = f'{current_app.config["tempdir"]}{filename}_constraints.csv'
    if os.path.exists(file):
        os.remove(file)
    else:
        print("The file does not exist")


# --- save


def save_poi_exceptions(filename, exceptions=[]):
    if len(exceptions) > 0:
        exceptions = pd.DataFrame(exceptions)
    else:
        exceptions = pd.DataFrame(columns=["x_col", "y_col", "x_coord", "y_coord", "radius"])
    exceptions.to_csv(f'{current_app.config["tempdir"]}{filename}_exceptions.csv', index=False)


def save_poi_constraints(filename, constraints=[{"col": cycle_column, "operator": "BETWEEN", "val1": 0, "val2": 100}]):
    constraints = pd.DataFrame(constraints)
    constraints.to_csv(f'{current_app.config["tempdir"]}{filename}_constraints.csv', index=False)


# --- load
def load_poi_exceptions(filename):
    path = f'{current_app.config["tempdir"]}{filename}_exceptions.csv'
    if not os.path.exists(path):  # if there is no constraints file yet, create a default one
        save_poi_exceptions(filename)

    df_exceptions = pd.read_csv(path)
    return df_exceptions


def load_poi_constraints(filename) -> pd.DataFrame:
    path = f'{current_app.config["tempdir"]}{filename}_constraints.csv'
    if not os.path.exists(path):  # if there is no constraints file yet, create a default one
        save_poi_constraints(filename)

    df_constraints = pd.read_csv(path)
    return df_constraints


# --- reset
@reaction_cime_api.route("/reset_poi_constraints/<filename>", methods=["GET"])
def reset_poi_constraints(filename):
    save_poi_constraints(filename)
    save_poi_exceptions(filename)  # TODO: do we want to reset exceptions too?
    return {"msg": "ok"}


# --- update
@reaction_cime_api.route("/add_poi_exceptions", methods=["OPTIONS", "POST"])
def add_poi_exceptions():
    if request.method == "POST":
        filename = request.form.get("filename")
        new_exceptions = json.loads(request.form["exceptions"])
        new_exceptions_df = pd.DataFrame(new_exceptions)

        if (
            "x_col" in new_exceptions_df.columns
            and "y_col" in new_exceptions_df.columns
            and "x_coord" in new_exceptions_df.columns
            and "y_coord" in new_exceptions_df.columns
            and "radius" in new_exceptions_df.columns
        ):

            exceptions_df = load_poi_exceptions(filename)
            exceptions_df = exceptions_df.append(new_exceptions_df, ignore_index=True)

            poi_count = (
                get_cime_dbo()
                .get_filter_mask(filename, get_poi_constraints_filter(filename, load_poi_constraints(filename), exceptions_df))["mask"]
                .sum()
            )  # check if constraints are limited enough
            save_poi_exceptions(filename, exceptions_df)
            if poi_count > MAX_POINTS:
                return {"msg": "Too many experiments within the selected filter. Only a subset of the data will be shown."}
            return {"msg": "ok"}

        if len(new_exceptions_df) <= 0:
            # save_poi_exceptions(filename)
            return {"msg": "ok"}

        return {
            "error": "wrong format; constraints must be of form [{'x_col': columnname, 'y_col': columnname, 'x_coord': number, 'y_coord': number, 'radius': number}]"
        }
    else:
        return {}


@reaction_cime_api.route("/update_poi_exceptions", methods=["OPTIONS", "POST"])
def update_poi_exceptions():
    if request.method == "POST":
        filename = request.form.get("filename")
        exceptions = json.loads(request.form["exceptions"])
        exceptions_df = pd.DataFrame(exceptions)

        if (
            "x_col" in exceptions_df.columns
            and "y_col" in exceptions_df.columns
            and "x_coord" in exceptions_df.columns
            and "y_coord" in exceptions_df.columns
            and "radius" in exceptions_df.columns
        ):
            poi_count = (
                get_cime_dbo()
                .get_filter_mask(filename, get_poi_constraints_filter(filename, load_poi_constraints(filename), exceptions_df))["mask"]
                .sum()
            )  # check if constraints are limited enough
            save_poi_exceptions(filename, exceptions_df)
            if poi_count > MAX_POINTS:
                return {"msg": "Too many experiments within the selected filter. Only a subset of the data will be shown."}
            return {"msg": "ok"}

        if len(exceptions_df) <= 0:
            save_poi_exceptions(filename)
            return {"msg": "ok"}

        return {
            "error": "wrong format; constraints must be of form [{'x_col': columnname, 'y_col': columnname, 'x_coord': number, 'y_coord': number, 'radius': number}]"
        }
    else:
        return {}


@reaction_cime_api.route("/update_poi_constraints", methods=["OPTIONS", "POST"])
def update_poi_constraints():
    if request.method == "POST":
        filename = request.form.get("filename")
        constraints = json.loads(request.form["constraints"])
        constraints_df = pd.DataFrame(constraints)

        if (
            "col" in constraints_df.columns
            and "operator" in constraints_df.columns
            and "val1" in constraints_df.columns
            and "val2" in constraints_df.columns
        ):
            poi_count = (
                get_cime_dbo()
                .get_filter_mask(filename, get_poi_constraints_filter(filename, constraints_df, load_poi_exceptions(filename)))["mask"]
                .sum()
            )  # check if constraints are limited enough
            if poi_count <= 0:
                return {"msg": "No experiments are found with this filter. Please adjust the settings."}
            save_poi_constraints(filename, constraints_df)
            if poi_count > MAX_POINTS:
                return {"msg": "Too many experiments within the selected filter. Only a subset of the data will be shown."}
            return {"msg": "ok"}

        if len(constraints_df) <= 0:
            save_poi_constraints(filename, pd.DataFrame(columns=["col", "operator", "val1", "val2"]))
            return {"msg": "ok"}

        return {
            "error": "wrong format; constraints must be of form [{'col': columnname, 'operator': 'BETWEEN'|'EQUALS', 'val1': firstValue, 'val2': secondValue}]"
        }
    else:
        return {}


# --- load and download
@reaction_cime_api.route("/get_poi_exceptions/<filename>", methods=["GET"])
def get_poi_exceptions(filename):
    exceptions_df = load_poi_exceptions(filename)
    return json.dumps(exceptions_df.to_dict("records")).encode("utf-8")


@reaction_cime_api.route("/get_poi_constraints/<filename>", methods=["GET"])
def get_poi_constraints(filename):
    constraint_df = load_poi_constraints(filename)
    return json.dumps(constraint_df.to_dict("records")).encode("utf-8")


@reaction_cime_api.route("/download_poi_exceptions/<filename>", methods=["GET"])
def download_poi_exceptions(filename):
    path = f'{current_app.config["tempdir"]}{filename}_exceptions.csv'
    return send_file(path, as_attachment=True)


@reaction_cime_api.route("/download_poi_constraints/<filename>", methods=["GET"])
def download_poi_constraints(filename):
    path = f'{current_app.config["tempdir"]}{filename}_constraints.csv'
    return send_file(path, as_attachment=True)


# --- filter
def map_constraint_operator(row):
    if row.operator == "BETWEEN":
        return f'"{row.col}" BETWEEN {row.val1} AND {row.val2}'
    if row.operator == "EQUALS":
        return f'"{row.col}"="{row.val1}"'
    return ""


def map_exceptions_filter(row):
    # TODO: make circle?
    x_lower = row.x_coord - row.radius
    x_upper = row.x_coord + row.radius
    y_lower = row.y_coord - row.radius
    y_upper = row.y_coord + row.radius
    return f'("{row.x_col}" BETWEEN {x_lower} AND {x_upper}) AND ("{row.y_col}" BETWEEN {y_lower} AND {y_upper})'


def get_poi_constraints_filter(filename, df_constraints: Optional[pd.DataFrame] = None, df_exceptions: Optional[pd.DataFrame] = None):
    if df_constraints is None:
        df_constraints = load_poi_constraints(filename)

    if df_exceptions is None:
        df_exceptions = load_poi_exceptions(filename)

    cols = list(set(df_constraints["col"]))
    col_filters = []
    for col in cols:
        df_constraints_col = df_constraints[df_constraints["col"] == col]

        # concatenate constraints of a single column with OR
        col_filter = " OR ".join(df_constraints_col.apply(map_constraint_operator, axis=1).tolist())
        col_filters.append(f"({col_filter})")

    constraints_filter_string = "true"  # if no filters are set, everything should be selected
    if len(col_filters) > 0:
        # concatenate constraints of different columns with AND
        constraints_filter_string = " AND ".join(col_filters)

    exceptions_filter_string = "false"  # if no exceptions are made, we default to false
    if len(df_exceptions) > 0:
        exceptions_filter_string = " OR ".join(df_exceptions.apply(map_exceptions_filter, axis=1).tolist())

    filter_string = f"({constraints_filter_string}) OR ({exceptions_filter_string})"

    print(filter_string)

    return filter_string


# --- retrieve


@reaction_cime_api.route("/get_poi_csv/<filename>", methods=["GET"])
def get_points_of_interest(filename):
    # domain = pd.read_csv("./temp-files/%s"%filename)
    # domain = get_cime_dbo().get_dataframe_from_table(filename)
    start_time = time.time()

    poi_domain, is_subsample = get_poi_df_from_db(filename, get_cime_dbo())  # TODO: tell front-end that it is subsampled

    if len(poi_domain) > 0:
        poi_domain = preprocess_dataset(poi_domain)

    csv_buffer = StringIO()
    poi_domain.to_csv(csv_buffer, index=False)

    delta_time = time.time() - start_time
    print("--- took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), "to load file %s" % filename)
    print("--- took %i min %f s to load file %s" % (delta_time / 60, delta_time % 60, filename))
    return csv_buffer.getvalue()


def get_poi_df_from_db(filename, cime_dbo):
    poi_domain, is_subsample = cime_dbo.get_dataframe_from_table_filter(
        filename, get_poi_constraints_filter(filename), max_datapoints=MAX_POINTS
    )
    return poi_domain, is_subsample


def get_poi_mask(filename, cime_dbo):
    mask = cime_dbo.get_filter_mask(filename, get_poi_constraints_filter(filename))
    return mask


# endregion

# region --------------- Parallel Coordinates ------------------


@reaction_cime_api.route("/get_csv_by_columns/<filename>", methods=["GET"])
def get_csv_by_columns(filename):
    columns = request.args.getlist("cols")
    domain = handle_dataset_cache(filename, columns)
    domain = domain[columns]

    csv_buffer = StringIO()
    domain.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


# endregion

# region --------------- nearest neighbors ------------------


@reaction_cime_api.route("/get_k_nearest_from_csv/<filename>/<x>/<y>/<k>", methods=["GET"])
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
    nearest_points_domain = get_cime_dbo().get_dataframe_from_table_complete_filter(
        filename, "ORDER BY ((x-(" + str(x) + "))*(x-(" + str(x) + ")))+((y-(" + str(y) + "))*(y-(" + str(y) + "))) LIMIT " + str(k) + ";"
    )

    response = make_response(nearest_points_domain.to_csv(index=False))
    response.headers["Content-Disposition"] = "attachment; filename=%s_kNN_%s_%s.csv" % (filename, x, y)
    response.headers["Content-type"] = "text/csv"
    return response


@reaction_cime_api.route("/get_radius_from_csv/<filename>/<x>/<y>/<d>", methods=["GET"])
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
    nearest_points_domain, _ = get_cime_dbo().get_dataframe_from_table_filter(
        filename, "((x-(" + str(x) + "))*(x-(" + str(x) + ")))+((y-(" + str(y) + "))*(y-(" + str(y) + "))) < " + str(r2)
    )

    response = make_response(nearest_points_domain.to_csv(index=False))
    response.headers["Content-Disposition"] = "attachment; filename=%s_kNN_%s_%s.csv" % (filename, x, y)
    response.headers["Content-type"] = "text/csv"
    return response


# endregion

# region --------------- aggregate dataset ------------------

# deprecated
@reaction_cime_api.route("/get_agg_csv/<filename>/<col_name>", methods=["GET"])
def get_aggregated_dataset(filename, col_name):
    range = {
        "x_min": request.args.get("x_min"),
        "x_max": request.args.get("x_max"),
        "y_min": request.args.get("y_min"),
        "y_max": request.args.get("y_max"),
    }
    filter = "x > {x_min} and x < {x_max} and y > {y_min} and y < {y_max}".format(**range)
    agg_domain, _ = get_cime_dbo().get_dataframe_from_table_filter(filename, filter, columns=["x", "y", col_name])
    agg_df = aggregate_by_col_interpolate(agg_domain, col_name, sample_size=200)  # TODO: dynamic sample_size

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


# deprecated
@reaction_cime_api.route("/get_agg_csv_cached/<filename>", methods=["GET"])
def get_aggregated_dataset_cached(filename):

    # value_col_name = request.args.get("value_col")
    # uncertainty_col_name = request.args.get("uncertainty_col")
    retrieve_cols = request.args.getlist("retrieve_cols")
    cache_cols = request.args.getlist(
        "cache_cols"
    )  # if value col or uncertainty col are not cached, we use this list of columns to prepare the new cached dataset
    sample_size = request.args.get("sample_size", default=200, type=int)
    range = {
        "x_min": float(request.args["x_min"]),
        "x_max": float(request.args["x_max"]),
        "y_min": float(request.args["y_min"]),
        "y_max": float(request.args["y_max"]),
    }

    agg_domain = handle_dataset_cache(filename, retrieve_cols + cache_cols)
    agg_domain = agg_domain[
        (agg_domain["x"] < range["x_max"])
        * (agg_domain["x"] > range["x_min"])
        * (agg_domain["y"] < range["y_max"])
        * (agg_domain["y"] > range["y_min"])
    ]

    agg_df = aggregate_by_col_interpolate(agg_domain, retrieve_cols, sample_size=sample_size)
    # agg_df = aggregate_by_col(agg_domain, retrieve_cols, sample_size=sample_size)

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


@reaction_cime_api.route("/get_hex_agg/<filename>/<x_channel>/<y_channel>", methods=["GET"])
def get_hexagonal_aggregation(filename, x_channel="x", y_channel="y"):

    retrieve_cols = request.args.getlist("retrieve_cols")
    aggregation_methods = request.args.getlist("aggregation_methods")
    cache_cols = request.args.getlist(
        "cache_cols"
    )  # if value col or uncertainty col are not cached, we use this list of columns to prepare the new cached dataset
    sample_size = 20  # request.args.get("sample_size", default=20, type=int)
    range = {
        "x_min": float(request.args["x_min"]),
        "x_max": float(request.args["x_max"]),
        "y_min": float(request.args["y_min"]),
        "y_max": float(request.args["y_max"]),
    }

    agg_domain = handle_dataset_cache(filename, retrieve_cols + cache_cols, x_channel=x_channel, y_channel=y_channel)
    agg_domain = agg_domain[
        (agg_domain[x_channel] < range["x_max"])
        * (agg_domain[x_channel] > range["x_min"])
        * (agg_domain[y_channel] < range["y_max"])
        * (agg_domain[y_channel] > range["y_min"])
    ]

    agg_df, wrong_points = hex_aggregate_by_col(
        agg_domain, retrieve_cols, aggregation_methods, range=None, sample_size=sample_size, x_channel=x_channel, y_channel=y_channel
    )  # TODO: what works better: range set to the boundaries of the dataset, or range set to the boundaries of the screen i.e. range=range

    wrong_df = wrong_points[list(set(retrieve_cols + [x_channel, y_channel]))]
    wrong_df["hex"] = False

    agg_df["hex"] = True

    agg_df = agg_df.append(wrong_df, ignore_index=True)

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


@reaction_cime_api.route("/get_value_range/<filename>/<col_name>", methods=["GET"])
def get_value_range(filename, col_name):
    range = get_cime_dbo().get_value_range_from_table(filename, col_name).iloc[0]
    return {"min": float(range["min"]), "max": float(range["max"])}


@reaction_cime_api.route("/get_category_values/<filename>/<col_name>", methods=["GET"])
def get_category_values(filename, col_name):
    df = handle_dataset_cache(filename, [col_name])
    distinct_vals = list(set(df[col_name]))
    return {"values": distinct_vals}


@reaction_cime_api.route("/get_category_count/<filename>/<col_name>", methods=["GET"])
def get_category_count(filename, col_name):
    # cat_count = get_cime_dbo().get_category_count(filename, col_name)
    # return json.dumps(cat_count.to_dict('records')).encode('utf-8')
    df = handle_dataset_cache(filename, [col_name])
    df = df[[col_name]]
    df["count"] = 0  # dummy column, such that it can be aggregated by count
    cat_count = df.groupby(col_name, as_index=False).count()
    return json.dumps(cat_count.to_dict("records")).encode("utf-8")


@reaction_cime_api.route("/get_category_count_of_hex/<filename>/<col_name>/<x_channel>/<y_channel>", methods=["GET"])
def get_category_count_of_hex(filename, col_name, x_channel, y_channel):
    hex_x = float(request.args["x"])
    hex_y = float(request.args["y"])
    hex_circ_radius = float(request.args["circ_radius"])

    # df = get_cime_dbo().get_dataframe_from_table(filename, columns=[x_channel, y_channel, col_name])
    df = handle_dataset_cache(filename, [x_channel, y_channel, col_name])

    window = is_inside_hex(
        df[x_channel], df[y_channel], hex_x=hex_x, hex_y=hex_y, radius=circ_radius_to_radius(hex_circ_radius), circ_radius=hex_circ_radius
    )
    df = df[window][[col_name]]
    df["count"] = 0  # dummy column, such that it can be aggregated by count
    cat_count = df.groupby(col_name, as_index=False).count()
    return json.dumps(cat_count.to_dict("records")).encode("utf-8")


@reaction_cime_api.route("/get_density/<filename>/<col_name>", methods=["GET"])
def get_density(filename, col_name):
    # data = get_cime_dbo().get_dataframe_from_table(filename, columns=[col_name])[col_name]
    data = handle_dataset_cache(filename, [col_name])[col_name]

    from scipy.stats import gaussian_kde

    density = gaussian_kde(data)
    data_min = min(data)
    data_max = max(data)
    x_vals = np.linspace(data_min, data_max, 100)
    y_vals = density(x_vals)

    # norm_sd like calculated in PSE: CoralDetail.tsx > getNormalizedSTD
    data_range = data_max - data_min
    if data_range == 0:
        data_range = 1
    return {"x_vals": list(x_vals), "y_vals": list(y_vals), "norm_sd": np.std((data - data_min) / data_range)}


@reaction_cime_api.route("/get_density_of_hex/<filename>/<col_name>/<x_channel>/<y_channel>", methods=["GET"])
def get_density_of_hex(filename, col_name, x_channel="x", y_channel="y"):
    hex_x = float(request.args["x"])
    hex_y = float(request.args["y"])
    hex_circ_radius = float(request.args["circ_radius"])

    # df = get_cime_dbo().get_dataframe_from_table(filename, columns=[x_channel, y_channel, col_name])
    df = handle_dataset_cache(filename, [x_channel, y_channel, col_name])

    window = is_inside_hex(
        df[x_channel], df[y_channel], hex_x=hex_x, hex_y=hex_y, radius=circ_radius_to_radius(hex_circ_radius), circ_radius=hex_circ_radius
    )
    data = df[window][col_name]

    from scipy.stats import gaussian_kde

    density = gaussian_kde(data)
    data_min = min(data)
    data_max = max(data)
    x_vals = np.linspace(data_min, data_max, 100)
    y_vals = density(x_vals)

    # norm_sd like calculated in PSE: CoralDetail.tsx > getNormalizedSTD
    data_range = data_max - data_min
    if data_range == 0:
        data_range = 1
    return {"x_vals": list(x_vals), "y_vals": list(y_vals), "norm_sd": np.std((data - data_min) / data_range)}


# endregion

# region --------------- projection ------------------


# depricated
@reaction_cime_api.route("/project_dataset", methods=["OPTIONS", "POST"])
def project_dataset():
    filename = request.form.get("filename")
    params = json.loads(request.form["params"])
    selected_feature_info = json.loads(request.form["selected_feature_info"])
    proj_df = get_cime_dbo().get_dataframe_from_table(filename, columns=list(selected_feature_info.keys()))

    # TODO: if params["useSelection"]: only project selected points --> do this in front-end? do this at all?

    # handle custom initialization of coordinates
    initialization = None
    if params["seeded"]:  # take custom seed from current coordinates
        initialization = get_poi_df_from_db(filename, get_cime_dbo())[0][["x", "y"]].values

    # rescale numerical values and encode categorical values
    proj_df, categorical_features = rescale_and_encode(proj_df, params, selected_feature_info)

    # handle custom metrics
    metric = params["distanceMetric"]
    if metric is None or metric == "":
        metric = "euclidean"
    if metric == "gower":  # we precompute the similarity matrix with the gower metric
        metric = "precomputed"
        normalized_values = gower.gower_matrix(proj_df, cat_features=categorical_features)
    else:  # otherwise, the similarity can be done by a pre-configured function and we just hand over the values
        normalized_values = proj_df.values

    # project the data
    if params["embedding_method"] == "umap":
        if initialization is None:
            initialization = "spectral"

        proj = umap.UMAP(
            n_neighbors=int(params["nNeighbors"]),
            n_components=2,
            metric=metric,
            n_epochs=int(params["iterations"]),
            init=initialization,
            verbose=True,
        )  # output_metric="euclidean" --> ? learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option
        proj_data = proj.fit_transform(normalized_values)
    elif params["embedding_method"] == "tsne":
        if initialization is None:
            initialization = "pca"
        proj = TSNE(
            2,
            metric=metric,
            perplexity=int(params["perplexity"]),
            n_iter=int(params["iterations"]),
            initialization=initialization,
            verbose=True,
        )  # , learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option
        proj_data = proj.fit(normalized_values)
    else:
        print("--- defaulting to PCA embedding")
        pca = PCA(n_components=2)
        proj_data = pca.fit_transform(normalized_values)

    # update the coordinates in the dataset
    start_time = time.time()
    get_cime_dbo().update_row_bulk(filename, proj_df.index, {"x": proj_data[:, 0], "y": proj_data[:, 1]})  # type: ignore
    delta_time = time.time() - start_time
    print("--- took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), "to update database")
    print("--- took %i min %f s to update database" % (delta_time / 60, delta_time % 60))

    print("return poi positions")
    df = get_poi_df_from_db(filename, get_cime_dbo())[0][["x", "y"]].to_dict("records")
    return {"done": "true", "embedding": df}


@reaction_cime_api.route("/terminate_projection_thread/<filename>", methods=["GET"])
def terminate_projection_thread(filename):
    current_app.config["TERMINATE_PROJECTION_" + filename] = True
    return {"msg": "ok"}


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
                proj_data, proj_df = self.run_rm_overlap()
            else:
                metric, normalized_values, proj_df, initialization = self.preprocess_dataset()

                self.msg = "calc embedding..."
                # project the data
                # --- umap
                if self.params["embedding_method"] == "umap":
                    proj_data = self.run_umap(metric, normalized_values, initialization)

                # --- tsne
                elif self.params["embedding_method"] == "tsne":
                    proj_data = self.run_tsne(metric, normalized_values, initialization)

                # --- pca
                else:
                    print("--- defaulting to PCA embedding")
                    pca = PCA(n_components=2)
                    proj_data = pca.fit_transform(normalized_values)

            self.current_step = self.params["iterations"]
            self.msg = "save embedding..."
            # update the coordinates in the dataset
            start_time = time.time()
            self.cime_dbo.update_row_bulk(self.filename, proj_df.index, {"x": proj_data[:, 0], "y": proj_data[:, 1]})  # type: ignore
            reset_dataset_cache(self.filename)  # reset cached data
            delta_time = time.time() - start_time
            print("--- took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), "to update database")
            print("--- took %i min %f s to update database" % (delta_time / 60, delta_time % 60))

            self.msg = "finished!"

        except Exception as e:
            self.msg = "error during projection"
            print("--Exception:--")
            print(e)

        finally:
            self.done = True

    def preprocess_dataset(self):
        self.msg = "load dataset..."
        proj_df = self.cime_dbo.get_dataframe_from_table(self.filename, columns=list(self.selected_feature_info.keys()))

        self.msg = "seeding..."
        # handle custom initialization of coordinates
        initialization = None
        if self.params["seeded"]:  # take custom seed from current coordinates
            # initialization = get_poi_df_from_db(self.filename, self.cime_dbo)[["x","y"]].values # does this make sense?
            initialization = proj_df[["x", "y"]].values

        # TODO: if params["useSelection"]: only project selected points --> do this in front-end? do this at all?

        self.msg = "rescale and encode..."
        # rescale numerical values and encode categorical values
        proj_df, categorical_features = rescale_and_encode(proj_df, self.params, self.selected_feature_info)

        self.msg = "calc metric..."
        # handle custom metrics
        metric = self.params["distanceMetric"]
        if metric is None or metric == "":
            metric = "euclidean"

        if metric == "gower":  # we precompute the similarity matrix with the gower metric
            metric = "precomputed"
            normalized_values = gower.gower_matrix(proj_df, cat_features=categorical_features)

            if initialization == "pca" and self.params["embedding_method"] == "tsne":
                print("--- calculating PCA for tSNE with Gowers distance")
                pca = PCA(n_components=2)
                initialization = pca.fit_transform(proj_df.values)

        else:  # otherwise, the similarity can be done by a pre-configured function and we just hand over the values
            normalized_values = proj_df.values

        return metric, normalized_values, proj_df, initialization

    def run_tsne(self, metric, normalized_values, initialization):
        # https://github.com/pavlin-policar/openTSNE/blob/master/openTSNE/callbacks.py
        from openTSNE.callbacks import Callback

        class CustomCallback(Callback):
            def __init__(self, parent, poi_mask):
                self.parent = parent
                self.poi_mask = poi_mask

            def optimization_about_to_start(self):
                """This is called at the beginning of the optimization procedure."""
                print("start optimization")
                pass

            def __call__(self, iteration, error, embedding):  # error is current KL divergence of the given embedding
                # poi_mask = get_poi_mask(self.parent.filename, self.parent.cime_dbo)
                self.parent.current_step = iteration
                if self.poi_mask is not None:
                    self.parent.emb = [{"x": row[0], "y": row[1]} for row in embedding[self.poi_mask]]

        #         return True # with this optimization will be interrupted

        if initialization is None:
            initialization = "pca"

        # instead of "None" give the poi_mask as parameter, if you want point updates in front-end: get_poi_mask(self.filename, self.cime_dbo)["mask"]
        # however, there is an error: Working outside of application context. (maybe due to changes in "get_poi_mask"?)
        proj = TSNE(
            2,
            metric=metric,
            callbacks_every_iters=10,
            callbacks=[CustomCallback(self, None)],
            perplexity=int(self.params["perplexity"]),
            n_iter=int(self.params["iterations"]),
            initialization=initialization,
            verbose=True,
        )  # , learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option
        proj_data = proj.fit(normalized_values)
        return proj_data

    def run_umap(self, metric, normalized_values, initialization):
        if initialization is None:
            initialization = "spectral"

        # https://pypi.org/project/tqdm/#parameters
        # use TQDM to update steps for front-end
        class CustomTQDM:
            def __init__(self, parent):
                self.parent = parent

            def write(self, str):
                try:
                    int(str)  # if it is not possible to parse it as int, we reached the end?
                    self.parent.current_step = str
                except ValueError:
                    print(str)

            def flush(self):
                pass

        proj = umap.UMAP(
            n_neighbors=int(self.params["nNeighbors"]),
            n_components=2,
            metric=metric,
            n_epochs=int(self.params["iterations"]),
            init=initialization,
            tqdm_kwds={"bar_format": "{n_fmt}", "file": CustomTQDM(self)},
            verbose=True,
        )  # output_metric="euclidean" --> ? learning_rate=float(params["learningRate"]) --> if we have user input for this, we would need "auto" as an option
        proj_data = proj.fit_transform(normalized_values)
        return proj_data

    def run_rm_overlap(self):
        self.msg = "load dataset..."
        proj_df = self.cime_dbo.get_dataframe_from_table(self.filename, columns=["x", "y"])

        self.msg = "calc embedding..."
        from .dgrid import DGrid

        # TODO: make parameters dynamic
        icon_width = 1
        icon_height = 1
        delta = 6
        start_time = time.time()

        from .dgrid_callback import CustomCallback as DGridCallback

        class MyCallback(DGridCallback):
            def __init__(self, parent):
                self.parent = parent

            def __call__(self, msg):
                # TODO: self.parent.current_step
                self.parent.msg = msg

        proj_data = DGrid(icon_width=icon_width, icon_height=icon_height, delta=delta, callbacks=[MyCallback(self)]).fit_transform(
            proj_df[["x", "y"]].values
        )
        delta_time = time.time() - start_time
        print("--- took", time.strftime("%H:%M:%S", time.gmtime(delta_time)), "calculate remove overlap")
        return proj_data, proj_df

    def get_id(self):

        # returns id of the respective thread
        if hasattr(self, "_thread_id"):
            return self._thread_id  # type: ignore
        for id, thread in threading._active.items():  # type: ignore
            if thread is self:
                return id

    def raise_exception(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, ctypes.py_object(SystemExit))
        if res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
            print("Exception raise failure")

        self.done = True
        self.msg = "client-side termination"


# differnt approach could have been with sockets: https://www.shanelynn.ie/asynchronous-updates-to-a-webpage-with-flask-and-socket-io/
@reaction_cime_api.route("/project_dataset_async", methods=["OPTIONS", "POST"])
def project_dataset_async():
    filename = request.form["filename"]
    params = json.loads(request.form["params"])
    selected_feature_info = json.loads(request.form["selected_feature_info"])
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
            yield json.dumps({"step": proj.current_step, "msg": proj.msg, "emb": proj.emb}).encode("utf-8")
            if terminate:
                proj.raise_exception()
                proj.join()
            if proj.done:  # proj.current_step >= params["iterations"] or proj.interrupt:
                df = get_poi_df_from_db(filename, cime_dbo)[0][["x", "y"]].to_dict("records")
                yield json.dumps({"step": None, "msg": None, "emb": df}).encode("utf-8")
                break

    # https://flask.palletsprojects.com/en/1.1.x/patterns/streaming/
    return Response(stream_with_context(my_generator()))


# endregion

# region --------------- chem specific -------------------


@reaction_cime_api.route("/get_mol_img", methods=["OPTIONS", "POST"])
def smiles_to_img_post():
    if request.method == "POST":
        smiles = request.form.get("smiles")
        img = smiles_to_base64(smiles)
        return {"data": img}
    else:
        return {}


@reaction_cime_api.route("/get_common_mol_img", methods=["OPTIONS", "POST"])
def smiles_list_to_common_substructure_img():
    if request.method == "POST":
        smiles_list = request.form.getlist("smiles_list")
        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}
        if len(smiles_list) == 1:
            ret = smiles_to_base64(smiles_list[0])
            return {"data": ret, "smiles": smiles_list[0]}

        mol_lst = []
        error_smiles = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)  # type: ignore
            if mol:
                mol_lst.append(mol)
            else:
                error_smiles.append(smiles)

        m = get_mcs(mol_lst)
        pil_img = Draw.MolToImage(m)

        if pil_img:
            buffered = BytesIO()
            pil_img.save(buffered, format="JPEG")  # type: ignore
            img_str = base64.b64encode(buffered.getvalue())
            return {"data": img_str.decode("utf-8"), "smiles": Chem.MolToSmiles(m)}  # type: ignore
        return {}
    else:
        return {}


@reaction_cime_api.route("/get_substructure_count", methods=["OPTIONS", "POST"])
def smiles_list_to_substructure_count():
    if request.method == "POST":
        smiles_list = request.form["smiles_list"].split(",")
        filter_smiles = request.form.get("filter_smiles")

        if len(smiles_list) == 0:
            return {"error": "empty SMILES list"}

        patt = Chem.MolFromSmiles(filter_smiles)  # type: ignore
        if patt:
            substructure_counts = [
                (smiles, len(Chem.MolFromSmiles(smiles).GetSubstructMatch(patt)))  # type: ignore
                for smiles in smiles_list
                if Chem.MolFromSmiles(smiles) is not None  # type: ignore
            ]
            return {"substructure_counts": substructure_counts}
        return {"error": "invalid SMILES filter"}
    else:
        return {}


# endregion

# region --------- clustering ---------


@reaction_cime_api.route("/segmentation", methods=["OPTIONS", "POST"])
def segmentation():
    if request.method == "POST":
        # clusterVal = request.forms.get("clusterVal")
        min_cluster_size_arg = request.form.get("min_cluster_size")
        min_cluster_samples_arg = request.form.get("min_cluster_samples")
        allow_single_cluster_arg = request.form.get("allow_single_cluster")
        x = request.form.get("X")
        assert x is not None
        x = np.array(x.split(","), dtype=np.float64)[:, np.newaxis].reshape((-1, 2))

        # many small clusters
        min_cluster_size = 5
        min_cluster_samples = 1
        allow_single_cluster = False

        if min_cluster_size_arg:
            min_cluster_size = int(min_cluster_size_arg)
        if min_cluster_samples_arg:
            min_cluster_samples = int(min_cluster_samples_arg)
        if allow_single_cluster_arg == "true":
            allow_single_cluster = bool(allow_single_cluster_arg)

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_cluster_samples,
            # prediction_data=True, # needed for soft clustering, or if we want to add points to the clustering afterwards
            allow_single_cluster=allow_single_cluster,  # maybe disable again
        )

        clusterer.fit_predict(x)

        # print(clusterer.labels_)
        # clusterer.probabilities_ = np.array(len(x))

        return {"result": [int(label) for label in clusterer.labels_]}

    else:
        return {}


# endregion
