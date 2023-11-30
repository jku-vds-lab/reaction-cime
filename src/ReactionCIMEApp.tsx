import * as React from 'react';
import {
  PSEContextProvider,
  API,
  Application,
  PluginRegistry,
  setItemLabel,
  setStoryLabel,
  setStoryBookLabel,
  setStoryTellingLabel,
  usePSESelector,
  RootActions,
  useCancellablePromise,
} from 'projection-space-explorer';
import { connect, ConnectedProps, useDispatch } from 'react-redux';
import HelpIcon from '@mui/icons-material/Help';
import { IconButton, Tooltip } from '@mui/material';
import { useVisynAppContext, VisynApp, VisynHeader } from 'visyn_core/app';
import { Anchor, useMantineTheme } from '@mantine/core';

import { BrowserRouter, useSearchParams } from 'react-router-dom';
import { IClientConfig } from 'visyn_core/base';
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
import vdsLogo from './assets/jku-vds-lab-logo.svg';
import bayerLogo from './assets/bayer_logo.svg';
import { BackendCSVLoader } from './Overrides/Dataset/BackendCSVLoader';
import { TabDocumentation } from './Utility/TabDocumentation';
import { BuildInfoContent, BuildInfoLogos } from './Utility/HeaderCustomization';
import { DummyEmbeddingController } from './Overrides/Embeddings/DummyEmbeddingController';

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

const ApplicationWrapper = connector(({ setMouseMoveFn, setMouseClickFn, resetViews }: Props) => {
  const { clientConfig } = useVisynAppContext();
  const loadedDataset = usePSESelector((state) => state.dataset);
  const [searchParams, setSearchParams] = useSearchParams();
  const dispatch = useDispatch();
  const { cancellablePromise } = useCancellablePromise();
  const abortController = React.useMemo(() => new AbortController(), []);
  const hasStartedRef = React.useRef(false);

  const startProjection = (msg: string) => {
    if (msg === 'init') {
      resetViews();
    }
  };

  const embeddings = [
    {
      id: 'umapRemote',
      name: 'UMAP',
      settings: { nneighbors: true },
      embController: new RemoteEmbeddingController('umap', startProjection),
      description:
        'Performs Uniform Manifold Approximation (UMAP) on the whole dataset using the chosen feature columns. This projects the high-dimensional dataset to a two-dimensional space that will then be shown as a scatterplot.',
      tooltip:
        'Performs Uniform Manifold Approximation (UMAP) on the whole dataset using the chosen feature columns. This method scales better than t-SNE with an increasing number of points.',
    },
    {
      id: 'tsneRemote',
      name: 't-SNE',
      settings: { perplexity: true },
      embController: new RemoteEmbeddingController('tsne', startProjection),
      description:
        'Performs t-distributed stochastic neighbor embedding (t-SNE) on the whole dataset using the chosen feature columns. This projects the high-dimensional dataset to a two-dimensional space that will then be shown as a scatterplot.',
      tooltip: 'Performs t-distributed stochastic neighbor embedding (t-SNE) on the whole dataset using the chosen feature columns.',
    },
    {
      id: 'pcaRemote',
      name: 'PCA',
      settings: {},
      embController: new RemoteEmbeddingController('pca', startProjection),
      description:
        'Performs principal component analysis (PCA) on the whole dataset using the chosen feature columns. This linearly projects the high-dimensional dataset to a two-dimensional space that will then be shown as a scatterplot.',
      tooltip: 'Performs principal component analysis (PCA) on the whole dataset using the chosen feature columns.',
    },
    {
      id: 'rmOverlap',
      name: 'Overlap removal',
      settings: { hideSettings: true },
      embController: new RemoteEmbeddingController('rmOverlap', startProjection),
      tooltip: 'Removes overlapping points by moving them to the nearest non-overlapping position. This helps to reduce visual clutter.',
      description:
        'Removes overlapping points of the current projection by moving them to the nearest non-overlapping position. This is particularly useful for large datasets after a projection like t-SNE or UMAP has been triggered to reduce visual clutter.',
    },
  ];

  const onLoadProject = React.useCallback(
    (tableId: string) => {
      new BackendCSVLoader().resolvePath(
        { path: tableId, uploaded: true },
        (dataset) => {
          dispatch(RootActions.loadDataset(dataset));
        },
        cancellablePromise,
        null,
        abortController,
      );
    },
    [dispatch, abortController, cancellablePromise],
  );

  React.useEffect(() => {
    if (hasStartedRef.current) {
      return;
    }
    hasStartedRef.current = true;

    const tableId = searchParams.get('project');

    // If the application started with a tableId, try to load that
    if (tableId) {
      onLoadProject(tableId);
    }
  }, [searchParams, onLoadProject]);

  React.useEffect(() => {
    const urlId = loadedDataset?.info?.path;

    if (!hasStartedRef.current) {
      return;
    }

    setSearchParams((prev) => {
      if (urlId) {
        return new URLSearchParams({
          project: urlId,
        });
      }

      return prev;
    });
  }, [loadedDataset?.info, setSearchParams]);

  return (
    <Application
      config={{
        preselect: {
          initOnMount: false,
        },
        baseUrl: ReactionCIMEBackendFromEnv.baseUrl,
      }}
      features={{
        embeddings: clientConfig?.publicVersion
          ? [
              {
                id: 'noProjection',
                name: 'Dummy Projection',
                settings: { hideSettings: true },
                embController: new DummyEmbeddingController(),
                description:
                  'Performs a projection on the whole dataset using the chosen feature columns. This projects the high-dimensional dataset to a two-dimensional space that will then be shown as a scatterplot. Projections are disabled in the Demo version.',
                tooltip: 'Performs a Projection on the whole dataset using the chosen feature columns.',
              },
            ]
          : embeddings,
        showVisibleProjections: false,
        showTrailSettings: false,
        enableFeatureWeighing: true,
        detailViewSplitRatio: [60, 40],
      }}
      overrideComponents={{
        tabContainerPrefix: TabDocumentation,
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
        contextMenuItems: [AddRegionExceptionMenuItem, SetFiltersToItemFeatures],
        detailViews: [
          {
            name: 'LineUp',
            view: <LineUpContext key="lineup" />,
            settings: LineUpTabPanel,
            alwaysRender: true,
          },
          {
            name: 'Parallel Coordinates',
            view: <PacoContext key="paco" />,
            settings: PacoTabPanel,
            alwaysRender: true,
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
            description: 'Aggregated dataset that should be shown in the background',
            icon: ReactionCIMEIcons.Aggregate,
          },
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
  const context = React.useMemo(() => new API<AppState>(null, createCIMERootReducer()), []);

  React.useEffect(() => {
    context.store.dispatch(setItemLabel({ label: 'experiment', labelPlural: 'experiments' }));
    context.store.dispatch(setStoryLabel({ label: 'group sequence', labelPlural: 'group sequences' }));
    context.store.dispatch(setStoryBookLabel({ label: 'collection', labelPlural: 'collections' }));
    context.store.dispatch(setStoryTellingLabel({ label: 'group comparison' }));
  }, [context]);

  const { user } = useVisynAppContext();
  const { clientConfig } = useVisynAppContext();
  const theme = useMantineTheme();

  return (
    <PSEContextProvider context={context}>
      <VisynApp
        header={
          clientConfig?.publicVersion ? (
            // eslint-disable-next-line react/jsx-no-useless-fragment
            <></>
          ) : (
            <VisynHeader
              backgroundColor={theme.colors.dark[theme.fn.primaryShade()]}
              components={{
                aboutAppModal: {
                  content: <BuildInfoContent />,
                  customerLogo: <BuildInfoLogos />,
                },
                beforeRight: (
                  <>
                    <Tooltip title="Opens the documentation page for CIME4R">
                      <IconButton color="primary" size="small" href="https://github.com/jku-vds-lab/reaction-cime#documentation-cime4r" target="_blank">
                        <HelpIcon fontSize="inherit" />
                      </IconButton>
                    </Tooltip>
                    <Anchor
                      href="https://www.bayer.com/"
                      rel="noreferrer"
                      target="_blank"
                      sx={{
                        // Center the image
                        display: 'flex',
                        alignItems: 'center',
                      }}
                    >
                      <img src={bayerLogo} alt="Bayer logo" style={{ height: '32px', transform: 'scale(1.1)' }} />
                    </Anchor>
                    <Anchor
                      href="https://jku-vds-lab.at/"
                      rel="noreferrer"
                      target="_blank"
                      sx={{
                        // Center the image
                        display: 'flex',
                        alignItems: 'center',
                      }}
                    >
                      <img src={vdsLogo} alt="JKU VDS Lab logo" style={{ height: '24px' }} />
                    </Anchor>
                  </>
                ),
              }}
            />
          )
        }
      >
        {user ? (
          <BrowserRouter>
            <ApplicationWrapper />
          </BrowserRouter>
        ) : null}
      </VisynApp>
    </PSEContextProvider>
  );
}
