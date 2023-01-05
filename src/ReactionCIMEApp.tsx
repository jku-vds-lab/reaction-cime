import * as React from 'react';
import { useEffect, useState } from 'react';
import { PSEContextProvider, API, Application, PluginRegistry, setItemLabel } from 'projection-space-explorer';
import { connect, ConnectedProps } from 'react-redux';
import { LineUpContext } from './LineUpContext';
import { LineUpTabPanel } from './Overrides/LineUpTabPanel';
import { AppState, CIME4RViewActions, createCIMERootReducer } from './State/Store';
import { AggregationTabPanel } from './Overrides/AggregationTabPanel';
import { DatasetTabPanel } from './Overrides/Dataset/DatasetTabPanel';
import { RemoteEmbeddingController } from './Overrides/Embeddings/RemoteEmbeddingController';
import { ReactionCIMEIcons } from './Utility/ReactionCIMEIcons';
import { HexAggregationLayer } from './Overrides/AggregationLayer/HexAggregationLayer';
import { setMouseClick, setMouseMove } from './State/MouseInteractionHooksDuck';
import { ReactionsPlugin } from './Overrides/Details/ReactionsPlugin';
import { ReactionCIMEBackendFromEnv } from './Backend/ReactionCIMEBackend';
import { PacoContext } from './PacoContext/PacoContext';
import { FilterTabPanel } from './Overrides/FilterTabPanel/FilterTabPanel';
import { PacoTabPanel } from './Overrides/PacoTabPanel/PacoTabPanel';
import { AddRegionExceptionMenuItem } from './Overrides/ContextMenu/AddRegionException';
import { SetFiltersToItemFeatures } from './Overrides/ContextMenu/SetFiltersToItemFeatures';
import 'react-data-grid/lib/styles.css';

PluginRegistry.getInstance().registerPlugin(new ReactionsPlugin());

const mapStateToProps = (state: AppState) => ({
  dataset_path: state.dataset?.info?.path,
  legendAttributes: state.genericFingerprintAttributes,
  globalLabels: state.globalLabels,
});

const mapDispatchToProps = (dispatch) => ({
  setMouseMoveFn: (value) => dispatch(setMouseMove(value)),
  setMouseClickFn: (value) => dispatch(setMouseClick(value)),
  resetViews: () => dispatch(CIME4RViewActions.resetViews()),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

const ApplicationWrapper = connector(({ setMouseMoveFn, dataset_path, setMouseClickFn, legendAttributes, globalLabels, resetViews }: Props) => {
  const startProjection = (msg: string) => {
    if (msg === 'init') {
      resetViews();
    }
  };

  useEffect(() => {
    if (legendAttributes != null) {
      // update cache in backend, when legendAttributes are changed
      const cols = legendAttributes.filter((item) => item.show).map((item) => item.feature);
      if (cols.length > 0) {
        ReactionCIMEBackendFromEnv.updateBackendCache(dataset_path, cols);
      }
    }
  }, [legendAttributes, dataset_path]);

  return (
    <Application
      config={{
        preselect: {
          initOnMount: false, // should default dataset be loaded? could specify url to default // TODO: define a default dataset that is already uploaded (e.g. domain.csv)
          // url: DATASETCONFIG[0].path
        },
        baseUrl: ReactionCIMEBackendFromEnv.baseUrl,
      }}
      features={{
        embeddings: [
          // {id:"umap", name:"UMAP", settings: DEFAULT_UMAP_SETTINGS},
          { id: 'umapRemote', name: 'UMAP', settings: { nneighbors: true }, embController: new RemoteEmbeddingController('umap', startProjection) },
          { id: 'tsneRemote', name: 't-SNE', settings: { perplexity: true }, embController: new RemoteEmbeddingController('tsne', startProjection) },
          { id: 'pcaRemote', name: 'PCA', settings: {}, embController: new RemoteEmbeddingController('pca', startProjection) },
          {
            id: 'rmOverlap',
            name: 'Overlap Removal',
            settings: { hideSettings: true },
            embController: new RemoteEmbeddingController('rmOverlap', startProjection),
          },
        ],
      }}
      overrideComponents={{
        mouseInteractionCallbacks: {
          onmousemove: (coords, event_used) => {
            setMouseMoveFn({ x: coords.x, y: coords.y, event_used });
          },
          onmouseclick: (coords, event_used, button) => {
            setMouseClickFn({ x: coords.x, y: coords.y, event_used, button });
          },
        },
        datasetTab: DatasetTabPanel,
        appBar: () => <div />,
        contextMenuItems: [
          // {key:"getkNN", title:"Download k-Nearest", function:(coords) => {
          //   handleBackgroundSelectionDownload(coords, dataset_path)
          // }},
          // {key:"addRegion", title:`Show ${globalLabels.itemLabelPlural} in this region`, function:(coords) => {
          //   // TODO
          //     ReactionCIMEBackendFromEnv.addPOIExceptions(dataset_path, [{x_col: "", y_col: "", x_coord: "", y_coord: "", radius: 10}]).then((res) => {
          //       if(res.msg === "ok"){
          //         // TODO: reload dataset
          //       }
          //     })
          // }}
          AddRegionExceptionMenuItem,
          SetFiltersToItemFeatures,
        ],
        detailViews: [
          {
            name: 'LineUp',
            view: <LineUpContext key="lineup" />,
            settings: LineUpTabPanel,
          },
          {
            name: 'Parallel Coordinates',
            view: <PacoContext key="paco" />,
            settings: PacoTabPanel,
          },
        ],
        tabs: [
          {
            name: 'filterDS',
            tab: FilterTabPanel,
            title: 'Filter',
            description: 'Filter subset of dataset that should be shown in the main visualization views',
            icon: ReactionCIMEIcons.Filter,
          },
          {
            name: 'aggregatDS',
            tab: AggregationTabPanel,
            title: 'Aggregate',
            description: 'Aggregated Dataset that should be shown in the background',
            icon: ReactionCIMEIcons.Aggregate,
          },
          // {
          //   name: "lineup",
          //   tab: <ViewsTabPanel overrideComponents={this.props.overrideComponents} splitRef={this.splitRef} />,
          //   title: "LineUp Integration",
          //   description: "Settings for LineUp Integration",
          //   icon: ReactionCIMEIcons.Table,
          // },
        ],
        layers: [
          {
            order: -1,
            component: HexAggregationLayer,
          },
        ],
      }}
    />
  );
});

export function ReactionCIMEApp() {
  // const [context] = useState(new API<RootState>(null, rootReducer))

  const [context] = useState(
    // new API<AppState>(null, createRootReducer(CIMEReducers))
    new API<AppState>(null, createCIMERootReducer()),
  );
  context.store.dispatch(setItemLabel({ label: 'experiment', label_plural: 'experiments' }));
  // context.store.dispatch(setDatasetEntriesAction(DATASETCONFIG))
  // context.store.getState().dataset...
  return (
    <PSEContextProvider context={context}>
      <ApplicationWrapper />
    </PSEContextProvider>
  );
}
