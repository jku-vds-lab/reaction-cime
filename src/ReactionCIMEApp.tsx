import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  createRootReducer,
  PSEIcons,
  DEFAULT_UMAP_SETTINGS,
} from "projection-space-explorer";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import { AggregationTabPanel } from "./Overrides/AggregationTabPanel";
import { AggregationLayer } from "./Overrides/AggregationLayer/AggregationLayer";
import { DatasetTabPanel } from "./Overrides/Dataset/DatasetTabPanel";
import { RemoteUMAPEmbeddingController } from "./Overrides/Embeddings/RemoteUMAPEmbeddingController";

export const DEMO = false;

// PluginRegistry.getInstance().registerPlugin(new ChemPlugin());

export function ReactionCIMEApp() {
  // const [context] = useState(new API<RootState>(null, rootReducer))

  const [context] = useState(
    new API<AppState>(null, createRootReducer(CIMEReducers))
  );

  // context.store.dispatch(setDatasetEntriesAction(DATASETCONFIG))
  

  return <PSEContextProvider context={context}><Application
    config={{
      preselect: {
        initOnMount: false, // should default dataset be loaded? could specify url to default // TODO: define a default dataset that is already uploaded (e.g. domain.csv)
        // url: DATASETCONFIG[0].path
      }
    }}
    features={{
      embeddings: [
        {id:"umap", name:"UMAP", settings: DEFAULT_UMAP_SETTINGS},
        // {id:"umapRemote", name:"UMAP Remote", embController: new RemoteUMAPEmbeddingController()}
      ],
    }}
    overrideComponents={{
      datasetTab: DatasetTabPanel,
      appBar: null,//CimeAppBar, --> remove when null
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
