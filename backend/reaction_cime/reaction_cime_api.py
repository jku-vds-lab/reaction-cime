import pickle
from flask import Blueprint, request
import logging
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from io import BytesIO
import base64

_log = logging.getLogger(__name__)

reaction_cime_api = Blueprint('reaction_cime', __name__)

@reaction_cime_api.route('/hello', methods=['GET'])
def get_uploaded_files_list():
    return "Hello World"

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

# ---helper functions---
def get_mcs(mol_list):

    if len(mol_list) <= 1:
        return Chem.MolFromSmiles("*")

    if type(mol_list[0]) == str:
        # TODO: handle invalid smiles
        mol_list = [Chem.MolFromSmiles(sm) for sm in mol_list]

    # completeRingsOnly=True # there are different settings possible here
    res = rdFMCS.FindMCS(mol_list, timeout=60, matchValences=False,
                         ringMatchesRingOnly=True, completeRingsOnly=True)
    if(res.canceled):
        patt = Chem.MolFromSmiles("*")
    else:
        patt = res.queryMol

    return patt


def smiles_to_base64(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m:
        return mol_to_base64(m)
    else:
        return "invalid smiles"


def mol_to_base64(m):
    pil_img = Draw.MolToImage(m)

    buffered = BytesIO()
    pil_img.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue())
    buffered.close()
    return img_str.decode("utf-8")
