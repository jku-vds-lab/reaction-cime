import base64
import re
from io import BytesIO

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS

# ---------------- preprocess dataset --------------------


def preprocess_dataset(domain):
    # calculates the error between measurement and prediction
    # target_column = [col for col in domain.columns if target_modifier in col.lower()][0]
    # if error_calc_modifier is not None:
    #     error_calc_col = [col for col in domain.columns if error_calc_modifier in col.lower()][0]
    #     error_col = domain.apply(lambda x: np.nan if x[cycle_column] < 0 else abs(x[target_column] - x['%s_mean_%i'%(error_calc_col, x[cycle_column])]), axis=1)
    # index = list(domain.columns).index(target_column)

    new_cols = generate_rename_list(domain)
    domain.columns = new_cols
    # if error_calc_modifier is not None:
    #     domain.insert(loc=index+1, column='error{"project":false,"paco":false,"real_column":false}', value=error_col)

    return domain


# def add_meta_info_time_series_data(column, timesteps, featureLabel, timeSeriesGroup, globalRange=None, colorMapping=None):
#     for t in range(timesteps):
#         modifiers = '"featureLabel":"%s", "timeSeriesGroup":"%s", "project":false'%(featureLabel, timeSeriesGroup)
#         if globalRange:
#             modifiers += ', "globalRange":%s'%globalRange
#         if colorMapping:
#             modifiers += ', "colorMapping":['
#             for c in colorMapping: # need to do it this way, because otherwise it gives single quotes...
#                 modifiers += '"%s",'%c
#             modifiers = modifiers[:-1] # remove last comma
#             modifiers += ']'

#         col = column%t
#         rename_dict[col] = '%s{%s}'%(col, modifiers)


def generate_global_ranges(domain):
    ranges = {}
    for col_sub in time_series_cols_diverging:
        global_min = min(domain[[col for col in domain.columns if col_sub in col]].fillna(0).min())
        global_max = max(domain[[col for col in domain.columns if col_sub in col]].fillna(0).max())

        global_range = max(abs(global_min), abs(global_max))
        # global_range = round(global_range) # do we want to be rounded?
        ranges[col_sub] = [-global_range, global_range]

    return ranges


def get_time_series_modifier(col, modifier, global_ranges):
    split_name = col.split("_")
    timestep = split_name.pop(-1)
    featureLabel = "_".join(split_name)

    # -- column is a tuple time series feature
    if len([ele for ele in time_series_tuples if ele in col]) > 0:
        featureVal = split_name.pop(-1)
        timeSeriesGroup = (
            "_".join(split_name) + ":" + featureVal
        )  # for pred_mean and pred_var we need to add a colon because lineup uses the colon to find time series with 2 variables
        modifier += '"paco":false,'  # show column in parallel coordinates
    else:
        timeSeriesGroup = "_".join(split_name)

    modifier += '"featureLabel":"%s", "timeSeriesGroup":"%s", "timestep":%s, "project":false' % (featureLabel, timeSeriesGroup, timestep)

    # -- column is a diverging time series feature
    ts_col_sub_lst = [ele for ele in time_series_cols_diverging if ele in col]
    if len(ts_col_sub_lst) > 0:
        ts_col_sub = ts_col_sub_lst[0]

        if global_ranges[ts_col_sub]:
            modifier += ', "globalRange":%s' % global_ranges[ts_col_sub]

        # add diverging colormap for diverging feature
        modifier += ', "colorMapping":['
        for c in ["#1e88e5", "#ffffff", "#ff0d57"]:  # need to do it this way, because otherwise it gives single quotes...
            modifier += '"%s",' % c
        modifier = modifier[:-1]  # remove last comma
        modifier += "]"

    return modifier


# cycle_column = "experimentCycle" # old version
cycle_column = "experiment_cycle"
# target_column = "yield" # old version
# target_column = "measured_yield"
target_modifier = "measured"
# error_calc_col = "pred" # old version
# error_calc_col = "predicted_yield"
# error_calc_modifier = "predicted"
time_series_tuples = ["predicted", "pred"]
time_series_cols_diverging = ["shap"]
smiles_modifier = "smiles"
# experiment_parameters = ["substrate_concentration", "sulfonyl_equiv", "base_equiv", "temperature"]
experiment_modifier = "exp_param"
# hide_lineup_summary_cols = ["sulfonyl_fluoride", "base", "solvent"]
hide_lineup_summary_modifier = "desc"

# error_calc_col = None
# target_column = "yield"
# cycle_column = "experimentCycle"
# time_series_tuples = ["pred"]
# time_series_cols_diverging = ["shap"]
# smiles_modifier = "smiles"
# experiment_parameters = ["concentration", "temperature", "Ligand_SMILES", "Base_SMILES", "Solvent_SMILES"]
# hide_lineup_summary_cols = [] #["reagent", "catalyst", "solvent"]


def generate_rename_list(domain):

    global_ranges = generate_global_ranges(domain)

    new_cols = []
    for col in domain.columns:
        modifier = ""
        col_name = col

        # -- column is a time series feature
        if re.search(r"_\d+$", col) is not None:
            modifier = get_time_series_modifier(col, modifier, global_ranges)
        elif experiment_modifier in col.lower():
            # elif col in experiment_parameters:
            modifier = '"project":true,"featureLabel":"exp_parameters","paco":true,"colName":"%s"' % col.replace(
                experiment_modifier, ""
            )  # TODO: use colName as visual name in frontend
        elif target_modifier in col.lower():
            # elif target_col == col:
            # col ends with _value bzw _step in lineup -> it belongs to a lineup time series
            # TODO: make this dynamic
            modifier = '"project":false,"paco":true,"lineup_meta_column":"pred_value"'  # signal lineup that it should add a meta_column with this label, that gives information for other columns
        elif col == cycle_column:
            # col ends with _value bzw _step in lineup -> it belongs to a lineup time series
            # TODO: make this dynamic
            modifier = '"project":false,"paco":false,"lineup_meta_column":"pred_step"'  # signal lineup that it should add a meta_column with this label, that gives information for other columns
        elif col == "groupLabel":
            groups = list(set(domain[domain[col] != "-1"][col]))
            groups.sort()
            clusterEdges = []
            indices = range(0, len(groups) - 1)
            for i in indices:
                clusterEdges.append([i, groups[i], groups[i + 1], ""])
            # need double quotes for javascript json parsing
            modifier = '"project":false,"paco":false,"edges":{"columns":["id","source","destination","name"],"index":%s,"data":%s}' % (
                str(list(indices)),
                str(clusterEdges).replace("'", '"'),
            )
        else:
            # -- column should not be shown in lineup or in the summary view
            # hide_col_list = [elem for elem in hide_lineup_summary_cols if elem in col]
            if hide_lineup_summary_modifier in col.lower():
                # if len(hide_col_list) > 0:
                modifier = '"noLineUp":true,"featureLabel":"%s","project":false,"colName":"%s"' % (
                    col.split("_")[0],
                    col.replace(hide_lineup_summary_modifier, ""),
                )

            # -- column does not have any special meaning
            else:
                modifier = '"project":false,"paco":false'

        # -- column is a smiles feature
        if smiles_modifier in col.lower():
            modifier += ',"imgSmiles":true'

        new_cols.append(
            '%s{"real_column":true,%s}' % (col_name, modifier)
        )  # "real_column" indicates that the column is actually in the dataset and was not derived

    return new_cols


# --- calculate aggregation of dataset


def get_grid_data(x, y, sample_size=200):
    # Create grid values
    xi = np.linspace(x.min(), x.max(), sample_size)
    yi = np.linspace(y.min(), y.max(), sample_size)
    Xi, Yi = np.meshgrid(xi, yi)

    return xi, yi, Xi, Yi


# deprecated
def aggregate_by_col_interpolate(df, value_cols, sample_size=200):
    from scipy.interpolate import griddata

    x = df["x"]
    y = df["y"]

    xi, yi, Xi, Yi = get_grid_data(x, y, sample_size)

    res_df = pd.DataFrame({"x": Xi.flatten(), "y": Yi.flatten()})

    # -----------------------
    # Interpolation on a grid
    # -----------------------
    for value_col in value_cols:
        z = df[value_col]
        zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method="linear")

        res_df[value_col] = zi.flatten()

    return res_df


# test if point is inside hexagon adapted from http://www.playchilla.com/how-to-check-if-a-point-is-inside-a-hexagon
def isInsideHex(points_x, points_y, hex_x, hex_y, radius, circ_radius):
    points_q2x = abs(points_x - hex_x)  # transform the test point locally and to quadrant 2
    points_q2y = abs(points_y - hex_y)  # transform the test point locally and to quadrant 2
    window = (points_q2x <= circ_radius) * (points_q2y <= radius)  # bounding test (since q2 is in quadrant 2 only 2 tests are needed)
    # if (points_q2x > circ_radius) or (points_q2y > radius):
    #     return False
    window *= (
        radius * circ_radius - radius * points_q2x - circ_radius / 2 * points_q2y
    ) >= 0  # finally the dot product can be reduced to this due to the hexagon symmetry
    return window


np_agg_methods_dict = {
    "min": np.nanmin,
    "max": np.nanmax,
    "mean": np.nanmean,
    "median": np.nanmedian,
}


def create_hex(df, hex_x, hex_y, radius, circ_radius, value_cols, aggregation_methods, xChannel, yChannel):
    window = isInsideHex(df[xChannel], df[yChannel], hex_x, hex_y, radius, circ_radius)
    window_df = df[window]
    if len(window_df) > 0:
        res = {xChannel: hex_x, yChannel: hex_y, "circ_radius": circ_radius}
        for i in range(len(value_cols)):
            value_col = value_cols[i]
            aggregation_method = aggregation_methods[i]

            res[value_col] = np_agg_methods_dict[aggregation_method](window_df[value_col])
        return res, window
    else:
        return None, window


def radius_to_circ_radius(radius):
    return 2 * radius / (3 ** (1 / 2))


def circ_radius_to_radius(circ_radius):
    return (3 ** (1 / 2)) * circ_radius / 2


def hex_aggregate_by_col(df, value_cols, aggregation_methods, range=None, sample_size=20, xChannel="x", yChannel="y"):
    x = df[xChannel]
    y = df[yChannel]

    if range is None:  # set range to be the maximum and minimum of the available dataset points
        x_min = x.min()
        x_max = x.max()
        y_min = y.min()
        y_max = y.max()
    else:  # set the range to be a custom range (e.g. if you zoom out, the number of hexagons should get smaller and smaller wrt screen space)
        x_min = range["x_min"]
        x_max = range["x_max"]
        y_min = range["y_min"]
        y_max = range["y_max"]

    x_delta = x_max - x_min
    y_delta = y_max - y_min
    delta = max(x_delta, y_delta)

    radius = (
        delta / sample_size / 2
    )  # distance from center to "flat" side -> r https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
    circ_radius = radius_to_circ_radius(
        radius
    )  # ((3*hex_radius*hex_radius)/4)**(1/2) # distance from center to corner (circumcircle radius); also length of one side in the hexagon -> a https://www-formula.com/geometry/radius-circumcircle/radius-circumcircle-regular-hexagon

    hexes = []

    test_points_used = np.zeros((len(df),))

    # make sure that each point is contained
    x_min -= circ_radius
    x_max += circ_radius
    y_min -= circ_radius
    y_max += circ_radius

    # even rows of hexes
    range_x = np.arange(x_min, x_max, circ_radius * 3)
    range_y = np.arange(y_min, y_max, radius * 2)
    for i in range_x:
        for j in range_y:
            hex, window = create_hex(df, i, j, radius, circ_radius, value_cols, aggregation_methods, xChannel, yChannel)
            test_points_used += window
            if hex is not None:
                hexes.append(hex)

    # odd rows of hexes
    range_x = np.arange(x_min - 1.5 * circ_radius, x_max + 1.5 * circ_radius, circ_radius * 3)
    range_y = np.arange(y_min - radius, y_max + radius, radius * 2)
    for i in range_x:
        for j in range_y:
            hex, window = create_hex(df, i, j, radius, circ_radius, value_cols, aggregation_methods, xChannel, yChannel)
            test_points_used += window
            if hex is not None:
                hexes.append(hex)

    # ---check that each point is used exactly once
    # print(test_points_used.min(), test_points_used.max()) # min and max have to be 1
    if (test_points_used - 1).sum() != 0:
        print(
            "--- attention! there is sth wrong with the point assignment of the hexagons (it might be that points are used several times, or points are not used at all) ---"
        )
    wrong_points = df[test_points_used != 1]  # TODO: remove debugging at some point
    return pd.DataFrame(hexes), wrong_points


# deprecated
def aggregate_by_col(df, value_cols, sample_size=20):
    x = df["x"]
    y = df["y"]
    xi, yi, Xi, Yi = get_grid_data(x, y, sample_size)  # (200,) (200,) (200, 200) (200, 200)
    # delta gives difference between two steps i.e. stepsize
    delta_x = (xi[1] - xi[0]) / 2
    delta_y = (yi[1] - yi[0]) / 2

    res_df = pd.DataFrame({"x": Xi.flatten(), "y": Yi.flatten()})

    for value_col in value_cols:
        z = df[value_col]

        zi = np.zeros((sample_size, sample_size))  # (200, 200)
        for i in range(sample_size):
            for j in range(sample_size):
                window = (
                    (df["x"] >= (xi[i] - delta_x))
                    * (df["x"] < (xi[i] + delta_x))
                    * (df["y"] >= (yi[j] - delta_y))
                    * (df["y"] < (yi[j] + delta_y))
                )
                window_vals = z[window]
                if len(window_vals) > 0:
                    zi[i, j] = np.nanmax(window_vals)
                else:
                    zi[i, j] = np.nan

        res_df[value_col] = zi.T.flatten()

    return res_df


# --- rescale and encode values
def rescale_and_encode(proj_df, params, selected_feature_info):
    categorical_feature_list = []
    for col in selected_feature_info.keys():
        info = selected_feature_info[col]

        if info["featureType"] == "String":
            proj_df = proj_df.drop(columns=[col])
            print("featureType: String --> TODO: handle")

        elif info["featureType"] == "Quantitative":
            if info["normalize"]:
                if params["normalizationMethod"] == "normalize01":  # scale values between [0;1]
                    # upper = info["range"]["max"] # do not use this! it is info from the front-end that only has POI dataset
                    # lower = info["range"]["min"] # do not use this! it is info from the front-end that only has POI dataset
                    upper = proj_df[col].max()
                    lower = proj_df[col].min()
                    div = upper - lower
                    if div == 0:  # when all values are equal in a column, the range is 0, which would lead to an error
                        div = 1
                    proj_df[col] = (proj_df[col] - lower) / div
                else:  # otherwise: "standardize" values to have 0 mean and unit standard deviation
                    mean = proj_df[col].mean()
                    std = proj_df[col].std()
                    if std <= 0:  # when all values are equal in a column, the standard deviation can be 0, which would lead to an error
                        std = 1
                    proj_df[col] = (proj_df[col] - mean) / std

        elif info["featureType"] == "Categorical":
            if params["encodingMethod"] == "onehot":
                hot_encoded = pd.get_dummies(proj_df[col], prefix=col, dummy_na=True)
                proj_df = proj_df.drop(columns=[col])
                proj_df = proj_df.join(hot_encoded)
            else:
                categorical_feature_list.append(col)
                proj_df[col] = pd.Categorical(proj_df[col]).codes

        elif info["featureType"] == "Date":
            proj_df = proj_df.drop(columns=[col])
            print("featureType: Date --> TODO: handle")

        elif info["featureType"] == "Binary":
            proj_df[col] = pd.Categorical(proj_df[col]).codes

        elif info["featureType"] == "Ordinal":
            proj_df = proj_df.drop(columns=[col])
            print("featureType: Ordinal --> TODO: handle")

        elif info["featureType"] == "Array":
            proj_df = proj_df.drop(columns=[col])
            print("featureType: Array --> TODO: handle")

    categorical_features = [col in categorical_feature_list for col in proj_df.columns]  # np.zeros((len(proj_df.columns),), dtype=np.bool)
    return proj_df, categorical_features


# ----------------- chem functions ----------------------
def get_mcs(mol_list):

    if len(mol_list) <= 1:
        return Chem.MolFromSmiles("*")

    if type(mol_list[0]) == str:
        # TODO: handle invalid smiles
        mol_list = [Chem.MolFromSmiles(sm) for sm in mol_list]

    # completeRingsOnly=True # there are different settings possible here
    res = rdFMCS.FindMCS(mol_list, timeout=60, matchValences=False, ringMatchesRingOnly=True, completeRingsOnly=True)
    if res.canceled:
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