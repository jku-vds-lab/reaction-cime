

# ---------------- preprocess dataset --------------------

import re

def add_meta_info_time_series_data(column, timesteps, featureLabel, lineUpGroup, globalRange=None, colorMapping=None):
    for t in range(timesteps):
        modifiers = '"featureLabel":"%s", "lineUpGroup":"%s", "project":false'%(featureLabel, lineUpGroup)
        if globalRange:
            modifiers += ', "globalRange":%s'%globalRange
        if colorMapping:
            modifiers += ', "colorMapping":['
            for c in colorMapping: # need to do it this way, because otherwise it gives single quotes...
                modifiers += '"%s",'%c
            modifiers = modifiers[:-1] # remove last comma
            modifiers += ']'
            
        col = column%t
        rename_dict[col] = '%s{%s}'%(col, modifiers)


def generate_global_ranges(domain):
    ranges = {}
    for col_sub in time_series_cols_diverging:
        global_min = min(domain[[col for col in domain.columns if col_sub in col]].min())
        global_max = max(domain[[col for col in domain.columns if col_sub in col]].max())

        global_range = max(abs(global_min), abs(global_max))
        # global_range = round(global_range) # do we want to be rounded?
        ranges[col_sub] = [-global_range,global_range]
        
    return ranges

def get_time_series_modifier(col, modifier, global_ranges):
    split_name = col.split("_")
    timestep = split_name.pop(-1)
    featureLabel = "_".join(split_name)

    # -- column is a tuple time series feature
    if len([ele for ele in time_series_tuples if ele in col]) > 0:
        lineUpGroup = ":".join(split_name) # for pred_mean and pred_var we need to add a colon because lineup uses the colon to find time series with 2 variables
    else: 
        lineUpGroup = "_".join(split_name)

    modifier += '"featureLabel":"%s", "lineUpGroup":"%s", "project":false'%(featureLabel, lineUpGroup)
    

    # -- column is a diverging time series feature
    ts_col_sub_lst = [ele for ele in time_series_cols_diverging if ele in col]
    if len(ts_col_sub_lst) > 0:
        ts_col_sub = ts_col_sub_lst[0]

        if global_ranges[ts_col_sub]:
            modifier += ', "globalRange":%s'%global_ranges[ts_col_sub]

        # add diverging colormap for diverging feature
        modifier += ', "colorMapping":['
        for c in ["#1e88e5","#ffffff","#ff0d57"]: # need to do it this way, because otherwise it gives single quotes...
            modifier += '"%s",'%c
        modifier = modifier[:-1] # remove last comma
        modifier += ']'

    return modifier


time_series_tuples = ["pred"]
time_series_cols_diverging = ["shap"]
smiles_modifier = "smiles"
experiment_parameters = ["substrate_concentration", "sulfonyl_equiv", "base_equiv", "temperature"]
hide_lineup_summary_cols = ["sulfonyl_fluoride", "base", "solvent"]

def generate_rename_list(domain):

    global_ranges = generate_global_ranges(domain)

    new_cols = []
    for col in domain.columns:
        modifier = ""
        col_name = col

        # -- column is a time series feature
        if re.search(r'_\d+$', col) != None:
            modifier = get_time_series_modifier(col, modifier, global_ranges)
        # -- column is a smiles feature
        elif smiles_modifier in col.lower():
            modifier = '"project":true,"imgSmiles":true,"featureLabel":"smiles"'
        elif col in experiment_parameters:
            modifier = '"project":true,"featureLabel":"exp_parameters"'
        else:
            # -- column should not be shown in lineup or in the summary view
            hide_col_list = [elem for elem in hide_lineup_summary_cols if elem in col]
            if len(hide_col_list) > 0:
                modifier = '"noLineUp":true,"featureLabel":"%s", "project":false'%hide_col_list[0]

            # -- column does not have any special meaning
            else:
                modifier = '"project":false'

        new_cols.append('%s{%s}'%(col_name, modifier))
    
    return new_cols


# --- calculate aggregation of dataset 
import pandas as pd
import numpy as np
def get_grid_data(x, y, sample_size=200):
    # Create grid values
    xi = np.linspace(x.min(), x.max(), sample_size)
    yi = np.linspace(y.min(), y.max(), sample_size)
    Xi, Yi = np.meshgrid(xi, yi)
    
    return xi, yi, Xi, Yi

def aggregate_col(df, value_col, sample_size=200):
    x = df["x"]
    y = df["y"]
    z = df[value_col]
    
    xi, yi, Xi, Yi = get_grid_data(x, y, sample_size)
    
    # -----------------------
    # Interpolation on a grid
    # -----------------------
    from scipy.interpolate import griddata
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')
    
    return pd.DataFrame({"x": Xi.flatten(), "y": Yi.flatten(), "val": zi.flatten()})



# ----------------- chem functions ----------------------
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from io import BytesIO
import base64

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
