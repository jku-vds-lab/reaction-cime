from typing import Any

import pandas as pd

from reaction_cime.helper_functions import rescale_and_encode


def test_rescale_and_encode():
    assert Any  # To avoid auto-removing the type
    df = pd.DataFrame(
        {
            # Numerical
            "num": [1, 2, 3, 4, 5],
            # Numerical
            "num_normalize": [1, 2, 3, 4, 5],
            # Categorical
            "cat": ["a", "b", "c", "d", "e"],
        }
    )  # type: Any

    params = {
        "perplexity": 50,
        "learningRate": 50,
        "nNeighbors": 15,
        "iterations": 500,
        "seeded": False,
        "useSelection": False,
        "method": "umapRemote",
        "distanceMetric": "euclidean",
        "normalizationMethod": "standardize",
        "encodingMethod": "onehot",
        "embedding_method": "umap",
    }
    features = {
        "num": {"featureType": "Quantitative", "weight": 1},
        "num_normalize": {"featureType": "Quantitative", "weight": 1, "normalize": True},
        "cat": {"featureType": "Categorical", "weight": 1},
    }

    featurizers = {}
    # For the numerical features, assume it stays the same
    categorical_features, weights, featurizer = rescale_and_encode(df, params, "num", features["num"])
    featurizers["num"] = featurizer
    featurizer(df, "num")
    assert (df["num"].values == [1, 2, 3, 4, 5]).all()

    # For the normalized numerical features, assume it got normalized
    categorical_features, weights, featurizer = rescale_and_encode(df, params, "num_normalize", features["num_normalize"])
    featurizers["num_normalize"] = featurizer
    featurizer(df, "num_normalize")
    normalized_num_values = df["num_normalize"].values
    assert normalized_num_values[2] == 0
    assert normalized_num_values[0] <= -1
    assert normalized_num_values[4] >= 1

    # For categorical features, assume it got one-hot encoded
    categorical_features, weights, featurizer = rescale_and_encode(df, params, "cat", features["cat"])
    featurizers["cat"] = featurizer
    featurizer(df, "cat")
    assert (df["cat_a"].values == [1, 0, 0, 0, 0]).all()
    assert (df["cat_b"].values == [0, 1, 0, 0, 0]).all()
    assert (df["cat_c"].values == [0, 0, 1, 0, 0]).all()
    assert (df["cat_d"].values == [0, 0, 0, 1, 0]).all()
    assert (df["cat_e"].values == [0, 0, 0, 0, 1]).all()

    # Now, test the case where we featurize a chunk of the original data and it should give the same results
    df = pd.DataFrame(
        {
            "num": [2, 3, 4],
            "num_normalize": [2, 3, 4],
            "cat": ["b", "c", "d"],
        }
    )

    featurizers["num"](df, "num")
    assert (df["num"].values == [2, 3, 4]).all()

    featurizers["num_normalize"](df, "num_normalize")
    assert (df["num_normalize"].values == normalized_num_values[1:4]).all()

    featurizers["cat"](df, "cat")
    assert (df["cat_a"].values == [0, 0, 0]).all()
    assert (df["cat_b"].values == [1, 0, 0]).all()
    assert (df["cat_c"].values == [0, 1, 0]).all()
    assert (df["cat_d"].values == [0, 0, 1]).all()
    assert (df["cat_e"].values == [0, 0, 0]).all()
