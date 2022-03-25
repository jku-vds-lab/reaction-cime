import { useEffect, useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  createRootReducer,
  PSEIcons,
  PluginRegistry,
} from "projection-space-explorer";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import { AggregationTabPanel } from "./Overrides/AggregationTabPanel";
import { DatasetTabPanel } from "./Overrides/Dataset/DatasetTabPanel";
import { RemoteEmbeddingController } from "./Overrides/Embeddings/RemoteEmbeddingController";
import { ReactionCIMEIcons } from "./Utility/ReactionCIMEIcons";
import { HexAggregationLayer } from "./Overrides/AggregationLayer/HexAggregationLayer";
import { handleBackgroundSelectionDownload } from "./Utility/Utils";
import { setMouseClick, setMouseMove } from "./State/MouseInteractionHooksDuck";
import { connect, ConnectedProps } from "react-redux";
import { ReactionsPlugin } from "./Overrides/Details/ReactionsPlugin";
import { ReactionCIMEBackendFromEnv } from "./Backend/ReactionCIMEBackend";
import { PacoContext } from "./PacoContext/PacoContext";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ReactionsPlugin());


export const ReactionCIMEApp = () => {
  // const [context] = useState(new API<RootState>(null, rootReducer))

  const [context] = useState(
    new API<AppState>(null, createRootReducer(CIMEReducers))
  );
  // context.store.dispatch(setDatasetEntriesAction(DATASETCONFIG))
  // context.store.getState().dataset...
  // context.store.dispatch(setMouseMove(coords))...
  
  return <PSEContextProvider context={context}><ApplicationWrapper></ApplicationWrapper></PSEContextProvider>

}

const mapStateToProps = (state: AppState) => ({
  dataset_path: state.dataset?.info?.path,
  legendAttributes: state.genericFingerprintAttributes,
});

const mapDispatchToProps = (dispatch) => ({
  setMouseMoveFn: (value) => dispatch(setMouseMove(value)),
  setMouseClickFn: (value) => dispatch(setMouseClick(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
};

const ApplicationWrapper = connector(({ setMouseMoveFn, dataset_path, setMouseClickFn, legendAttributes }: Props) => {
  useEffect(() => {
    if(legendAttributes != null){
      // update cache in backend, when legendAttributes are changed
      const cols = legendAttributes.filter((item) => item.show).map((item) => item.feature)
      if(cols.length > 0){
        ReactionCIMEBackendFromEnv.updateBackendCache(dataset_path, cols)
      }
    }
  }, [legendAttributes])

  return <Application
    config={{
      preselect: {
        initOnMount: false, // should default dataset be loaded? could specify url to default // TODO: define a default dataset that is already uploaded (e.g. domain.csv)
        // url: DATASETCONFIG[0].path
      }
    }}
    features={{
      embeddings: [
        // {id:"umap", name:"UMAP", settings: DEFAULT_UMAP_SETTINGS},
        {id:"umapRemote", name:"UMAP", settings: {nneighbors:true}, embController: new RemoteEmbeddingController("umap")},
        {id:"tsneRemote", name:"t-SNE", settings: {perplexity:true}, embController: new RemoteEmbeddingController("tsne")},
        {id:"pcaRemote", name:"PCA", settings: {}, embController: new RemoteEmbeddingController("pca")},
        {id:"rmOverlap", name:"Overlap Removal", settings: {hideSettings:true}, embController: new RemoteEmbeddingController("rmOverlap")}, //TODO: implement overplot removal in backend; tell PSE to not show any settings
      ],
    }}
    overrideComponents={{
      mouseInteractionHooks: {
        "mousemove": (coords, event_used) => {setMouseMoveFn({x: coords.x, y: coords.y, event_used: event_used})},
        // "mousedown": (coords) => {setMouseDownFn({x: coords.x, y: coords.y})},
        // "mouseup": (coords) => {setMouseUpFn({x: coords.x, y: coords.y})},
        "mouseclick": (coords, event_used, button) => { setMouseClickFn({x: coords.x, y: coords.y, event_used: event_used, button: button})} 
      },
      datasetTab: DatasetTabPanel,
      appBar: () => <div></div>,
      contextMenuItems: [{key:"getkNN", title:"Download k-Nearest", function:(coords) => {
        handleBackgroundSelectionDownload(coords, dataset_path)
      }}],
      detailViews: [ //TODO: switch between the two
        {
          name: "paco",
          // @ts-ignore
          view: <PacoContext key={"paco"}></PacoContext>
        },
        {
          name: "lineup",
          // @ts-ignore
          view: <LineUpContext key={"lineup"}></LineUpContext>
        },
      ],
      tabs: [
        {
          name: "aggregatDS",
          //@ts-ignore
          tab: AggregationTabPanel,
          title: "Aggregate",
          description: "Aggregated Dataset that should be shown in the background",
          icon: ReactionCIMEIcons.Aggregate
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
          component: () => <HexAggregationLayer></HexAggregationLayer>//<HexAggregationLayer></HexAggregationLayer>//<AggregationLayer></AggregationLayer> //<AggregationContourLayer></AggregationContourLayer>//
        }
       ]
    }}
  />
})

