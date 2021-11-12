import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  createRootReducer,
  PSEIcons,
  setDatasetEntriesAction,
} from "projection-space-explorer";
import * as THREE from 'three';
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import { DATASETCONFIG } from "./datasetconfig";
import { AggregationTabPanel } from "./Overrides/AggregationTabPanel";
import { GLHeatmap } from "./Overrides/AggregationLayer/GLHeatmap";
import { AggregationLayer } from "./Overrides/AggregationLayer/AggregationLayer";


// PluginRegistry.getInstance().registerPlugin(new ChemPlugin());

export function ReactionCIMEApp() {
  // const [context] = useState(new API<RootState>(null, rootReducer))

  const [context] = useState(
    new API<AppState>(null, createRootReducer(CIMEReducers))
  );

  context.store.dispatch(setDatasetEntriesAction(DATASETCONFIG))
  

  return <PSEContextProvider context={context}><Application
    config={{
      preselect: {
        initOnMount: true, // should default dataset be loaded? could specify url to default
        url: DATASETCONFIG[0].path
      }
    }}
    features={{
      disableEmbeddings: {
        tsne: true,
        forceatlas: true,
      },
    }}
    overrideComponents={{
      appBar: CimeAppBar,
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
      layers: [{
        order: -1,
        component: () => <AggregationLayer></AggregationLayer>
      }]
    }}
  /></PSEContextProvider>

}
