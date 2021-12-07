import { trackPromise } from "react-promise-tracker";
import {
  AVector,
  CSVLoader,
  Dataset,
  DatasetType,
  IVector,
  Loader,
  useCancellablePromise,
} from "projection-space-explorer";
import * as d3v5 from "d3v5";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";

function convertFromCSV(vectors) {
  return vectors.map((vector) => {
    return AVector.create(vector);
  });
}

export class BackendCSVLoader implements Loader {
  vectors: IVector[] = [];
  datasetType: DatasetType = DatasetType.None;

  loadingArea = "global_loading_indicator";

  resolvePath(entry: any, finished: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: AbortController) {
    if (entry.uploaded) {
      // TODO: Instead of localstorage store it in state?
      // localStorage.setItem("id", entry.path);
      // use file that is already uploaded to backend
      this.loadPOICSV(finished, entry, cancellablePromise, modifiers, controller);
    } else {
      trackPromise(
        fetch(entry.path, { signal: controller?.signal })
          .then((response) => response.blob())
          .then((result) =>
            this.resolveContent(
              result,
              finished,
              cancellablePromise,
              modifiers,
              controller
            )
          )
          .catch((error) => {
            console.log(error);
          }),
        this.loadingArea
      );
    }
  }

  resolveContent(file: any, onChange: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: AbortController) {
    const promise = cancellablePromise
      ? cancellablePromise(
        ReactionCIMEBackendFromEnv.upload_csv_file(file, controller), controller)
      : ReactionCIMEBackendFromEnv.upload_csv_file(file, controller);
    trackPromise(
      promise
        .then((uploaded) => {
          //console.log("UPLOADED", uploaded);
          this.loadPOICSV(onChange, { display: "", type: this.datasetType, path: uploaded.id }, cancellablePromise, modifiers, controller);
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loadingArea
    );
  }

  loadPOICSV(finished: (dataset: Dataset) => void, entry, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: AbortController) {
    // request the server to return a csv file using the unique filename
    const poi_path = ReactionCIMEBackendFromEnv.baseUrl + "/get_poi_csv/" + entry.path;
    const promise = cancellablePromise
      ? cancellablePromise(
          d3v5.csv(poi_path, {...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal,}),
          controller
        ) : d3v5.csv(poi_path, {...ReactionCIMEBackendFromEnv.fetchParams, signal: controller?.signal,});
    trackPromise(
      promise
        .then((vectors) => {
            if(vectors.length <= 0){
                console.log("'points of interest' dataset is empty");
                alert("'points of interest' dataset is empty")
            }else{
                this.vectors = convertFromCSV(vectors);
                this.datasetType = DatasetType.Chem;
                new CSVLoader().resolve(finished, this.vectors, this.datasetType, entry);
            }
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loadingArea
    );
  }

}
