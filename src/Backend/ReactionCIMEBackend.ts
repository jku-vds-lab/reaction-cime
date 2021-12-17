
import * as d3v5 from "d3v5";
import {useCancellablePromise} from "projection-space-explorer";
import { trackPromise } from "react-promise-tracker";

export class ReactionCIMEBackend {
  protected smiles_cache = {};
  protected smiles_highlight_cache = {};
  protected cache = {};

  constructor(
    public readonly baseUrl: string,
    public readonly fetchParams: Parameters<typeof fetch>[1] = {}
  ) {}

  protected handleSmilesCache = (smiles: string, highlight = false) => {
    //already downloaded this image -> saved in smiles cache
    if (highlight) {
      return this.smiles_highlight_cache[smiles];
    } else {
      return this.smiles_cache[smiles];
    }
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
    if (Object.keys(data).includes("error")) {
      alert(data["error"]);
    }
    return data;
  };



  public getUploadedFiles = async (): Promise<string[]> => {
    let path = this.baseUrl + "/get_uploaded_files_list";

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
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
    let path = this.baseUrl + "/delete_file/" + filename;

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert("file could not be deleted. please, try again");
        console.log(error);
      });
  };


  public getStructureFromSmiles = (
    id: string | number,
    smiles: string,
    highlight = false,
    controller
  ) => {
    const cached_data = this.handleSmilesCache(smiles, highlight);
    if (cached_data) {
      return this.async_cache(cached_data);
    }

    const formData = new FormData();
    formData.append("smiles", smiles);
    formData.append("filename", id?.toString());

    let path = this.baseUrl + "/get_mol_img";
    if (highlight) {
      path += "/highlight";
    }

    return fetch(path, {
      ...this.fetchParams,
      method: "POST",
      body: formData,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        this.setSmilesCache(smiles, highlight, data["data"]);
        return data["data"];
      })
      .catch((error) => {
        // alert("could not load structure");
        console.log(error);
      });
  };


  public getMCSFromSmilesList = async (formData: FormData, controller?) => {
    let my_fetch;
    if (controller) {
      my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
        method: "POST",
        body: formData,
        signal: controller?.signal,
      });
    } else {
      my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
        method: "POST",
        body: formData,
      });
    }
    return my_fetch
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((response) => response["data"])
      .catch((error) => {
        // alert("could not get maximum common substructure")
        console.log(error);
      });
  };

  public getSubstructureCount = async (smiles_list, filter) => {
    const formData = new FormData();
    formData.append("smiles_list", smiles_list);
    formData.append("filter_smiles", filter);
    return fetch(this.baseUrl + "/get_substructure_count", {
      method: "POST",
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .then((data) => {
        if (Object.keys(data).includes("substructure_counts"))
          return data["substructure_counts"];
        else throw Error("Backend responded with error: " + data["error"]);
      })
      .catch((error) => {
        alert("could not find substructure match");
        console.log(error);
      });
  };



  public calculateHDBScanClusters = async (
    X,
    min_cluster_size,
    min_cluster_samples,
    allow_single_cluster
  ) => {
    const formData = new FormData();
    formData.append("min_cluster_size", min_cluster_size);
    formData.append("min_cluster_samples", min_cluster_samples);
    formData.append("allow_single_cluster", allow_single_cluster);
    formData.append("X", X);
    return fetch(this.baseUrl + "/segmentation", {
      method: "POST",
      body: formData,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        alert("error when calculating clusters");
        console.log(error);
      });
  };


  public terminate_projection = async (filename): Promise<{ response: any }> => {
    console.log("backend: terminate_projection")
    let path = this.baseUrl + "/terminate_projection_thread/" + filename;

    return fetch(path, {
      ...this.fetchParams,
      method: "GET",
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        console.log(error);
      });
  };

  public project_dataset = async(filename:string, params:object, selected_feature_info:object, controller?) => {
    const formData = new FormData();
    formData.append("filename", filename);
    formData.append("params", JSON.stringify(params));
    formData.append("selected_feature_info", JSON.stringify(selected_feature_info));
    return fetch(this.baseUrl + "/project_dataset_async", {
      method: "POST",
      body: formData,
      signal: controller?.signal,
    }).then(this.handleErrors)
      .catch((error) => {
        console.log(error);
      })
  }


  public upload_csv_file = async (
    file,
    controller?
  ): Promise<{ name: string; id: number }> => {
    // upload the csv file to the server
    // the response is a unique filename that can be used to make further requests
    const formData_file = new FormData();
    formData_file.append("myFile", file);
    formData_file.append("file_size", file.size);
  
    const promise = fetch(this.baseUrl + "/upload_csv", {
      ...this.fetchParams,
      method: "POST",
      body: formData_file,
      signal: controller?.signal,
    })
      .then(this.handleErrors)
      .then((response) => response.json())
      .then(this.handleJSONErrors)
      .catch((error) => {
        if (error.name === "AbortError") {
          console.log("Fetch aborted");
        } else {
          alert("error when uploading file. it might be too big");
          console.log(error);
        }
      });
    return promise;
  };


  loadingArea = "global_loading_indicator";
  public loadAggCSV(finished: (dataset: any) => void, path, column, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: AbortController) {
    // request the server to return a csv file using the unique filename
    const agg_path = ReactionCIMEBackendFromEnv.baseUrl + "/get_agg_csv/" + path + "/" + column; // TODO: make dynamic
    const promise = cancellablePromise
      ? cancellablePromise(
          d3v5.csv(agg_path, {...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal,}), controller
        ) : d3v5.csv(agg_path, {...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal,});
    trackPromise(
      promise
        .then((vectors) => {
            if(vectors.length <= 0){
                console.log("aggregation dataset is empty");
                alert("aggregation dataset is empty")
            }else{
                finished(vectors)
            }
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loadingArea
    );
  }
}





// Use the environment variables defined in the .env file
if (!process.env.REACT_APP_CIME_BACKEND_URL) {
  console.error("The ENV-variable REACT_APP_CIME_BACKEND_URL must be set.");
}

export const ReactionCIMEBackendFromEnv = new ReactionCIMEBackend(
  process.env.REACT_APP_CIME_BACKEND_URL,
  {
    // credentials: "omit",
  }
);
