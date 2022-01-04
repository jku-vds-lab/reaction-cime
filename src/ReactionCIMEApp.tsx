import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  createRootReducer,
  PSEIcons,
} from "projection-space-explorer";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import { AggregationTabPanel } from "./Overrides/AggregationTabPanel";
import { AggregationLayer } from "./Overrides/AggregationLayer/AggregationLayer";
import { DatasetTabPanel } from "./Overrides/Dataset/DatasetTabPanel";
import { RemoteEmbeddingController } from "./Overrides/Embeddings/RemoteEmbeddingController";
import { ReactionCIMEBackendFromEnv } from "./Backend/ReactionCIMEBackend";

export const DEMO = false;

// PluginRegistry.getInstance().registerPlugin(new ChemPlugin());


export const ReactionCIMEApp = () => {
  // const [context] = useState(new API<RootState>(null, rootReducer))

  const [context] = useState(
    new API<AppState>(null, createRootReducer(CIMEReducers))
  );

  // context.store.dispatch(setDatasetEntriesAction(DATASETCONFIG))
  
//   <MenuItem onClick={() => {
//     var coords = CameraTransformations.screenToWorld(
//       {
//         x: this.mouseController.currentMousePosition.x,
//         y: this.mouseController.currentMousePosition.y,
//       },
//       this.createTransform()
//     );
//     console.log('Pressed "Download k-Nearest" option from context-menu with coords :>> ', coords);
//     this.props.setCimeBackgroundSelection(coords);
//     handleClose()
// }}>Download k-Nearest</MenuItem>

  return <PSEContextProvider context={context}><Application
    config={{
      preselect: {
        initOnMount: false, // should default dataset be loaded? could specify url to default // TODO: define a default dataset that is already uploaded (e.g. domain.csv)
        // url: DATASETCONFIG[0].path
      }
    }}
    features={{
      embeddings: [
        // {id:"umap", name:"UMAP", settings: DEFAULT_UMAP_SETTINGS},
        {id:"umapRemote", name:"UMAP Remote", settings: {nneighbors:true}, embController: new RemoteEmbeddingController("umap")},
        {id:"tsneRemote", name:"t-SNE Remote", settings: {perplexity:true}, embController: new RemoteEmbeddingController("tsne")},
        {id:"pcaRemote", name:"PCA Remote", settings: {}, embController: new RemoteEmbeddingController("pca")}
      ],
    }}
    overrideComponents={{
      datasetTab: DatasetTabPanel,
      appBar: null,//CimeAppBar, --> remove when null
      contextMenuItems: [{key:"getkNN", title:"Download k-Nearest", function:(coords) => {
        handleBackgroundSelectionDownload(coords, context.store.getState().dataset?.info?.path)
      }}],
      detailViews: [
        {
          name: "lineup",
          //@ts-ignore
          view: LineUpContext
        },
      ],
      tabs: [
        {
          name: "aggregatDS",
          //@ts-ignore
          tab: AggregationTabPanel,
          title: "Aggregate",
          description: "Aggregated Dataset that should be shown in the background",
          icon: PSEIcons.Dataset,
        },
        {
          name: "lineup",
          //@ts-ignore
          tab: LineUpTabPanel,
          title: "LineUp Integration",
          description: "Settings for LineUp Integration",
          icon: PSEIcons.PseLineup,
        },
      ],
      layers: [
        {
          order: -1,
          component: () => <AggregationLayer></AggregationLayer>
        }
       ]
    }}
  /></PSEContextProvider>

}

/**
 * This is merely a helper function to decompose the code into smaller individual segments.
 * It is used to handle the changing background selection prop, which might, e.g., be triggered when a user clicks the k-nearest neighbor option in the context menu.
 * Specifically, it checks whether the parameters are correct, and if so, sends a query to the db to fetch the k-nearest entires to the click, with k being defined by a textfield.
 * The response triggers the download of a csv with these entries.
 * Afterwards, the background selection is reset, to make sure other prop updates do not trigger this db query and download.
 * @param {any} coords - The prop that holds x and y coordinates of the clicks
 * @returns {void} - no return value
 */
 function handleBackgroundSelectionDownload(coords: any, filename: string) {
  // if input for checking k-nearest neighbors (x,y coordinates and k) are not undefined
  if (
    typeof coords?.x !== "undefined" &&
    typeof coords?.y !== "undefined" &&
    (document.getElementById("knn-textfield") as HTMLInputElement)?.value !==
      "undefined"
  ) {
    let k = +(document.getElementById("knn-textfield") as HTMLInputElement)
      ?.value;
    // if input k is neither integer nor below 1
    if (k < 1 || k % 1 !== 0) {
      // warn user
      alert("Invalid input for k-nearest neighbors.");
    } else {
      // otherwise send request to db and download response in browser
      ReactionCIMEBackendFromEnv.getkNearestData(
        filename,
        coords?.x,
        coords?.y,
        (document.getElementById("knn-textfield") as HTMLInputElement)?.value
      );
    }
  }
}
