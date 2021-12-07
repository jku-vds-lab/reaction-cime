import pickle
from flask import Blueprint, request, current_app, abort, jsonify
import logging
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from io import BytesIO
import base64
from .helper_functions import generate_rename_list, get_mcs, smiles_to_base64, aggregate_col

_log = logging.getLogger(__name__)

reaction_cime_api = Blueprint('reaction_cime', __name__)

@reaction_cime_api.route('/hello', methods=['GET'])
def hello():
    return "Hello World"


@reaction_cime_api.route('/get_uploaded_files_list', methods=['GET'])
def get_uploaded_files_list():
    return jsonify(list(get_cime_dbo().get_table_names()))

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

    # TODO: make dynamic, by which feature we want to filter (e.g. user could change the settings in the front-end maybe with parallel coordinates?)
    # poi_domain = domain[domain['yield']>0]
    print(filename)
    poi_domain = get_cime_dbo().get_dataframe_from_table_filter(filename, "yield > 0")
    print(len(poi_domain))

    new_cols = generate_rename_list(poi_domain)
    poi_domain.columns = new_cols

    csv_buffer = StringIO()
    poi_domain.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


@reaction_cime_api.route('/get_agg_csv/<filename>/<col_name>', methods=['GET'])
def get_aggregated_dataset(filename, col_name):
    agg_domain = get_cime_dbo().get_dataframe_from_table(filename, columns=["x", "y", col_name])
    agg_df = aggregate_col(agg_domain, col_name, sample_size=200) # TODO: dynamic sample_size

    csv_buffer = StringIO()
    agg_df.to_csv(csv_buffer, index=False)

    return csv_buffer.getvalue()


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

