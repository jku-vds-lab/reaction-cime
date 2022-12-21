import * as d3v5 from 'd3v5';
import { useCancellablePromise } from 'projection-space-explorer';
import { trackPromise } from 'react-promise-tracker';

export class ReactionCIMEBackend {
  protected smiles_cache = {};

  protected smiles_highlight_cache = {};

  protected cache = {};

  constructor(public readonly baseUrl: string, public readonly fetchParams: Parameters<typeof fetch>[1] = {}) {}

  protected handleSmilesCache = (smiles: string, highlight = false) => {
    // already downloaded this image -> saved in smiles cache
    if (highlight) {
      return this.smiles_highlight_cache[smiles];
    }
    return this.smiles_cache[smiles];
  };

  protected setSmilesCache = (smiles, highlight = false, data) => {
    if (highlight) this.smiles_highlight_cache[smiles] = data;
    else this.smiles_cache[smiles] = data;
  };

  protected async_cache = async (cached_data) => {
    return cached_data;
  };

  handleCache = (key) => {
    if (this.cache[key]) return Object.assign(this.cache[key]); // return copy of cached object
    return null;
  };

  setCache = (key, value) => {
    this.cache[key] = value;
  };

  handleErrors = (response) => {
    if (!response.ok) {
      throw Error(response.statusText);
    }
    return response;
  };

  handleJSONErrors = (data) => {
    if (Object.keys(data).includes('error')) {
      alert(data.error);
    }
    return data;
  };

  public getUploadedFiles = async (): Promise<string[]> => {
    const path = `${this.baseUrl}/get_uploaded_files_list`;

    return fetch(path, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        // alert("could not load uploaded filenames.")
        console.log(error);
      });
  };

  public deleteFile = async (filename): Promise<{ deleted: any }> => {
    const path = `${this.baseUrl}/delete_file/${encodeURIComponent(filename)}`;

    return fetch(path, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert('file could not be deleted. please, try again');
        console.log(error);
      });
  };

  public getStructureFromSmiles = (smiles: string, highlight = false, controller) => {
    const cached_data = this.handleSmilesCache(smiles, highlight);
    if (cached_data) {
      return this.async_cache(cached_data);
    }

    const formData = new FormData();
    formData.append('smiles', smiles);

    let path = `${this.baseUrl}/get_mol_img`;
    if (highlight) {
      path += '/highlight';
    }

    return fetch(path, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        this.setSmilesCache(smiles, highlight, data.data);
        return data.data;
      })
      .catch((error) => {
        // alert("could not load structure");
        console.log(error);
      });
  };

  public getMCSFromSmilesList = async (formData: FormData, controller?) => {
    let my_fetch;
    if (controller) {
      my_fetch = fetch(`${this.baseUrl}/get_common_mol_img`, {
        ...this.fetchParams,
        method: 'POST',
        body: formData,
        signal: controller?.signal,
      });
    } else {
      my_fetch = fetch(`${this.baseUrl}/get_common_mol_img`, {
        ...this.fetchParams,
        method: 'POST',
        body: formData,
      });
    }
    return my_fetch
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((response) => response.data)
      .catch((error) => {
        // alert("could not get maximum common substructure")
        console.log(error);
      });
  };

  public getSubstructureCount = async (smiles_list, filter) => {
    const formData = new FormData();
    formData.append('smiles_list', smiles_list);
    formData.append('filter_smiles', filter);
    return fetch(`${this.baseUrl}/get_substructure_count`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        if (Object.keys(data).includes('substructure_counts')) return data.substructure_counts;
        throw Error(`Backend responded with error: ${data.error}`);
      })
      .catch((error) => {
        alert('could not find substructure match');
        console.log(error);
      });
  };

  public calculateHDBScanClusters = async (X, min_cluster_size, min_cluster_samples, allow_single_cluster) => {
    const formData = new FormData();
    formData.append('min_cluster_size', min_cluster_size);
    formData.append('min_cluster_samples', min_cluster_samples);
    formData.append('allow_single_cluster', allow_single_cluster);
    formData.append('X', X);
    return fetch(`${this.baseUrl}/segmentation`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert('error when calculating clusters');
        console.log(error);
      });
  };

  public terminate_projection = async (filename): Promise<{ response: any }> => {
    const path = `${this.baseUrl}/terminate_projection_thread/${encodeURIComponent(filename)}`;

    return fetch(path, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  /**
   * Utilizes d3v5.csv() and the corresponding API route to issue an HTTP GET request for a CSV.
   * This CSV will be returned by a database query in which the k nearest entries to the coordinates (x,y) are selected.
   * @param {string} filename - The filename which serves as a link to the dataset, i.e., the table name in the database.
   * @param {string} x - The x coordinate of the point for which the nearest neighbours should be returned
   * @param {string} y - The y coordinate of the point for which the nearest neighbours should be returned
   * @param {string} k - The amount of neighbours to return
   * @returns {Promise} - Returns a promise of the vectors returned by the database query.
   */
  public getkNearestData = async (filename: string, x: string, y: string, k: string) => {
    // TODO consider accepting integer parameters and converting them here to append them to the path string
    // console.log('calling get_nearest_data: filename, x, y, k :>> ', filename, x, y, k);
    const path = `${this.baseUrl}/get_k_nearest_from_csv/${encodeURIComponent(filename)}/${x}/${y}/${k}`;
    // window.location.href = path;
    window.open(path);
  };

  public project_dataset = async (filename: string, params: object, selected_feature_info: object, controller?) => {
    const formData = new FormData();
    formData.append('filename', filename);
    formData.append('params', JSON.stringify(params));
    formData.append('selected_feature_info', JSON.stringify(selected_feature_info));
    return fetch(`${this.baseUrl}/project_dataset_async`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => {
        this.resetAggregationCache();
        return response;
      })
      .catch((error) => {
        console.log(error);
      });
  };

  public upload_csv_file = async (file, controller?): Promise<{ name: string; id: number }> => {
    // upload the csv file to the server
    // the response is a unique filename that can be used to make further requests
    const formData_file = new FormData();
    formData_file.append('myFile', file);
    formData_file.append('file_size', file.size);

    const promise = fetch(`${this.baseUrl}/upload_csv`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData_file,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        if (error.name === 'AbortError') {
          console.log('Fetch aborted');
        } else {
          alert('error when uploading file. it might be too big');
          console.log(error);
        }
      });
    return promise;
  };

  public loadNODatapoints = async (filename: string) => {
    return fetch(`${this.baseUrl}/get_no_datapoints/${encodeURIComponent(filename)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadValueRange = async (filename: string, col_name: string) => {
    return fetch(`${this.baseUrl}/get_value_range/${encodeURIComponent(filename)}/${encodeURIComponent(col_name)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadCategoryValues = async (filename: string, col_name: string) => {
    return fetch(`${this.baseUrl}/get_category_values/${encodeURIComponent(filename)}/${encodeURIComponent(col_name)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public updateBackendCache = async (filename: string, cache_cols: string[]) => {
    let cache_cols_string = '';
    if (cache_cols != null) {
      for (const key in cache_cols) {
        const col = cache_cols[key];
        cache_cols_string += `&cache_cols=${col}`;
      }
    }
    return fetch(`${this.baseUrl}/update_cache/${encodeURIComponent(filename)}?dummy=1${cache_cols_string}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadCategoryCount = async (filename: string, col_name: string) => {
    return fetch(`${this.baseUrl}/get_category_count/${encodeURIComponent(filename)}/${encodeURIComponent(col_name)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadCategoryCountOfHex = async (filename: string, col_name: string, xChannel: string, yChannel: string, x: number, y: number, circ_radius: number) => {
    if (xChannel == null) {
      xChannel = 'x';
    }
    if (yChannel == null) {
      yChannel = 'y';
    }

    return fetch(
      `${this.baseUrl}/get_category_count_of_hex/${encodeURIComponent(filename)}/${encodeURIComponent(
        col_name,
      )}/${xChannel}/${yChannel}?x=${x}&y=${y}&circ_radius=${circ_radius}`,
      {
        ...this.fetchParams,
        method: 'GET',
      },
    )
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadDensity = async (filename: string, col_name: string) => {
    return fetch(`${this.baseUrl}/get_density/${encodeURIComponent(filename)}/${encodeURIComponent(col_name)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadDensityOfHex = async (filename: string, col_name: string, xChannel: string, yChannel: string, x: number, y: number, circ_radius: number) => {
    if (xChannel == null) {
      xChannel = 'x';
    }
    if (yChannel == null) {
      yChannel = 'y';
    }

    return fetch(
      `${this.baseUrl}/get_density_of_hex/${encodeURIComponent(filename)}/${encodeURIComponent(
        col_name,
      )}/${xChannel}/${yChannel}?x=${x}&y=${y}&circ_radius=${circ_radius}`,
      {
        ...this.fetchParams,
        method: 'GET',
      },
    )
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  protected agg_dataset_cache: [{ x: { min: number; max: number }; y: { min: number; max: number }; data: any }] = null;

  protected cur_agg_path = '';

  // protected cur_agg_value_col = "";
  // protected cur_agg_uncertainty_col = "";
  // reset when reprojecting, when dataset is changed, and when column is changed
  protected resetAggregationCache = () => {
    this.agg_dataset_cache = null;
  };

  protected cache_agg_data(range: { x: { min: number; max: number }; y: { min: number; max: number } }, vectors: any) {
    if (this.agg_dataset_cache == null) {
      this.agg_dataset_cache = [{ ...range, data: vectors }];
    } else {
      this.agg_dataset_cache.push({ ...range, data: vectors });
    }
  }

  protected is_within_boundaries(val1: number, val2: number) {
    const bound = Math.abs(val2 * 0.5);
    if (val1 < val2 + bound && val1 > val2 - bound) {
      return true;
    }
    return false;
  }

  protected handleAggregationCache = (
    path: string,
    value_col: string,
    uncertainty_col: string,
    range: { x: { min: number; max: number }; y: { min: number; max: number } },
    xChannel?: string,
    yChannel?: string,
  ) => {
    if (this.agg_dataset_cache == null || this.agg_dataset_cache.length <= 0) return null;

    xChannel = xChannel == null ? 'x' : xChannel;
    xChannel = yChannel == null ? 'y' : yChannel;
    if (
      this.cur_agg_path !== path ||
      !(value_col in this.agg_dataset_cache[0].data) ||
      !(uncertainty_col in this.agg_dataset_cache[0].data) ||
      !(xChannel in this.agg_dataset_cache[0].data) ||
      !(yChannel in this.agg_dataset_cache[0].data)
    ) {
      // if(this.cur_agg_path !== path || this.cur_agg_value_col !== value_col || this.cur_agg_uncertainty_col !== uncertainty_col){
      this.resetAggregationCache();
      this.cur_agg_path = path;
      // this.cur_agg_value_col = value_col;
      // this.cur_agg_uncertainty_col = uncertainty_col;
      return null;
    }

    const filtered = this.agg_dataset_cache.filter((value) => {
      // if it lies within a certain range of the cached data, we take it
      if (
        this.is_within_boundaries(range.x.min, value.x.min) &&
        this.is_within_boundaries(range.x.max, value.x.max) &&
        this.is_within_boundaries(range.y.min, value.y.min) &&
        this.is_within_boundaries(range.y.max, value.y.max)
      ) {
        return true;
      }
      return false;
    });
    if (filtered.length === 1) return filtered[0];
    if (filtered.length > 1) {
      // TODO: sort by similarity and return most similar
      // filtered.sort()
      return filtered[0];
    }

    return null;
  };

  public loadHexAgg(
    finished: (dataset: any) => void,
    path: string,
    xChannel: string,
    yChannel: string,
    value_column: string,
    uncertainty_col: string,
    cache_cols: string[],
    sample_size: number,
    aggregationMethod: any,
    range: { x: { min: number; max: number }; y: { min: number; max: number } },
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    controller?: AbortController,
    loadingArea?: string,
  ) {
    if (xChannel == null) {
      xChannel = 'x';
    }
    if (yChannel == null) {
      yChannel = 'y';
    }

    const cached_data = this.handleAggregationCache(path, value_column, uncertainty_col, range, xChannel, yChannel);
    let promise = null;
    if (cached_data) {
      console.log('cached');
      promise = this.async_cache(cached_data.data);
    } else {
      let retrieve_cols = `retrieve_cols=${value_column}`;
      retrieve_cols += `&aggregation_methods=${aggregationMethod.valueAggregationMethod}`;
      if (uncertainty_col !== 'None' && uncertainty_col != null) {
        retrieve_cols += `&retrieve_cols=${uncertainty_col}`;
        retrieve_cols += `&aggregation_methods=${aggregationMethod.uncertaintyAggregationMethod}`;
      }

      let cache_cols_string = '';
      if (cache_cols != null) {
        for (const key in cache_cols) {
          const col = cache_cols[key];
          cache_cols_string += `&cache_cols=${col}`;
        }
      }

      let sample_size_str = '';
      if (sample_size != null) {
        sample_size_str = `&sample_size=${sample_size}`;
      }

      if (range == null) {
        range = { x: { min: -10000, max: 10000 }, y: { min: -10000, max: 10000 } };
      }
      const range_string = `&x_min=${range.x.min}&x_max=${range.x.max}&y_min=${range.y.min}&y_max=${range.y.max}`;

      const agg_path = `${ReactionCIMEBackendFromEnv.baseUrl}/get_hex_agg/${path}/${xChannel}/${yChannel}?${retrieve_cols}${cache_cols_string}${sample_size_str}${range_string}`;

      promise = cancellablePromise
        ? cancellablePromise(d3v5.csv(agg_path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal }), controller)
        : d3v5.csv(agg_path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal });
    }
    trackPromise(
      promise
        .then((vectors) => {
          if (vectors.length <= 0) {
            console.log('aggregation dataset is empty');
            alert('aggregation dataset is empty');
          } else {
            if (cached_data == null) {
              this.cache_agg_data(range, vectors);
            }
            finished(vectors);
          }
        })
        .catch((error) => {
          console.log(error);
        }),
      loadingArea,
    );
  }

  public loadAggCSV(
    finished: (dataset: any) => void,
    path: string,
    value_column: string,
    uncertainty_col: string,
    cache_cols: string[],
    sample_size: number,
    range: { x: { min: number; max: number }; y: { min: number; max: number } },
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    controller?: AbortController,
    loadingArea?: string,
  ) {
    if (range == null) {
      range = { x: { min: -10000, max: 10000 }, y: { min: -10000, max: 10000 } };
    }

    const cached_data = this.handleAggregationCache(path, value_column, uncertainty_col, range);
    let promise = null;
    if (cached_data) {
      promise = this.async_cache(cached_data.data);
    } else {
      let retrieve_cols = `retrieve_cols=${value_column}`;
      if (uncertainty_col !== 'None' && uncertainty_col != null) {
        retrieve_cols += `&retrieve_cols=${uncertainty_col}`;
      }

      let cache_cols_string = '';
      if (cache_cols != null) {
        for (const key in cache_cols) {
          const col = cache_cols[key];
          cache_cols_string += `&cache_cols=${col}`;
        }
      }

      let sample_size_str = '';
      if (sample_size != null) {
        sample_size_str = `&sample_size=${sample_size}`;
      }

      const range_string = `&x_min=${range.x.min}&x_max=${range.x.max}&y_min=${range.y.min}&y_max=${range.y.max}`;

      // request the server to return a csv file using the unique filename
      // const agg_path = ReactionCIMEBackendFromEnv.baseUrl + "/get_agg_csv/" + path + "/" + column + "?x_min=" + range.x.min + "&x_max=" + range.x.max + "&y_min=" + range.y.min + "&y_max="+range.y.max; // TODO: make dynamic
      const agg_path = `${ReactionCIMEBackendFromEnv.baseUrl}/get_agg_csv_cached/${path}?${retrieve_cols}${cache_cols_string}${sample_size_str}${range_string}`; // TODO: make dynamic

      promise = cancellablePromise
        ? cancellablePromise(d3v5.csv(agg_path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal }), controller)
        : d3v5.csv(agg_path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal });
    }
    trackPromise(
      promise
        .then((vectors) => {
          if (vectors.length <= 0) {
            console.log('aggregation dataset is empty');
            alert('aggregation dataset is empty');
          } else {
            if (cached_data == null) {
              this.cache_agg_data(range, vectors);
            }
            finished(vectors);
          }
        })
        .catch((error) => {
          console.log(error);
        }),
      loadingArea,
    );
  }

  public loadPacoCSV(
    finished: (dataset: any) => void,
    filename: string,
    cols: string[],
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    controller?: AbortController,
    loadingArea?: string,
  ) {
    let cols_string = '?dummy=0';
    for (const i in cols) {
      const col = cols[i];
      cols_string += `&cols=${col}`;
    }

    const path = `${ReactionCIMEBackendFromEnv.baseUrl}/get_csv_by_columns/${encodeURIComponent(filename)}${cols_string}`;

    const promise = cancellablePromise
      ? cancellablePromise(d3v5.csv(path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal }), controller)
      : d3v5.csv(path, { ...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal });

    trackPromise(
      promise
        .then((vectors) => {
          if (vectors.length <= 0) {
            console.log('dataset is empty');
            alert('dataset is empty');
          } else {
            finished(vectors);
          }
        })
        .catch((error) => {
          console.log(error);
        }),
      loadingArea,
    );
  }

  public loadPOIExceptions = async (filename: string) => {
    return fetch(`${this.baseUrl}/get_poi_exceptions/${encodeURIComponent(filename)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public loadPOIConstraints = async (filename: string) => {
    return fetch(`${this.baseUrl}/get_poi_constraints/${encodeURIComponent(filename)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public resetPOIConstraints = async (filename: string) => {
    return fetch(`${this.baseUrl}/reset_poi_constraints/${encodeURIComponent(filename)}`, {
      ...this.fetchParams,
      method: 'GET',
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public addPOIExceptions = async (filename: string, exceptions: any[]) => {
    const formData = new FormData();
    formData.append('exceptions', JSON.stringify(exceptions));
    formData.append('filename', filename);
    return fetch(`${this.baseUrl}/add_poi_exceptions`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public updatePOIExceptions = async (filename: string, exceptions: any[]) => {
    const formData = new FormData();
    formData.append('exceptions', JSON.stringify(exceptions));
    formData.append('filename', filename);
    return fetch(`${this.baseUrl}/update_poi_exceptions`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public updatePOIConstraints = async (filename: string, constraints: any[]) => {
    const formData = new FormData();
    formData.append('constraints', JSON.stringify(constraints));
    formData.append('filename', filename);
    return fetch(`${this.baseUrl}/update_poi_constraints`, {
      ...this.fetchParams,
      method: 'POST',
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public downloadPOIConstraints = async (filename: string) => {
    const path = `${this.baseUrl}/download_poi_constraints/${encodeURIComponent(filename)}`;
    // window.location.href = path;
    window.open(path);
  };

  public uploadPOIConstraints = async (filename: string, file) => {
    const contents = await file.text();
    const data = d3v5.csvParse(contents);
    return ReactionCIMEBackendFromEnv.updatePOIConstraints(filename, data);
    // var reader = new FileReader();

    // reader.onload = function(e) {
    //   var contents = e.target.result;
    //   var data = d3v5.csvParse(contents);
    //   return ReactionCIMEBackendFromEnv.updatePOIConstraints(filename, data)
    // };

    // reader.readAsText(file);
  };
}

// Use the environment variables defined in the .env file
if (process.env.REACT_APP_CIME_BACKEND_URL == null) {
  console.error('The ENV-variable REACT_APP_CIME_BACKEND_URL must be set.');
}

export const ReactionCIMEBackendFromEnv = new ReactionCIMEBackend(process.env.REACT_APP_CIME_BACKEND_URL, {
  // credentials: "omit",
});
