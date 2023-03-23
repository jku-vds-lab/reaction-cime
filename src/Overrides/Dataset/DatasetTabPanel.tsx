import * as React from 'react';
import { Accordion, AccordionDetails, AccordionSummary, Box, Button, Divider, Grid, Tooltip, Typography } from '@mui/material';
import { ExpandMore, FileUpload, GetApp, InfoOutlined } from '@mui/icons-material';
import { Dataset, IProjection, isEntityId, RootActions, useCancellablePromise, UtilityActions } from 'projection-space-explorer';
import { connect, ConnectedProps, ReactReduxContext, useDispatch } from 'react-redux';
import { useState } from 'react';
import ManageSearchIcon from '@mui/icons-material/ManageSearch';
import { useVisynAppContext } from 'visyn_core';
import { DatasetDrop } from './DatasetDrop';
import { AppState, CIME4RViewActions } from '../../State/Store';
import { UploadedFiles } from './UploadedFiles';
import { BackendCSVLoader } from './BackendCSVLoader';
import { setTriggerUpdate as setTriggerUpdateAction } from '../../State/HandleDatasetDuck';
import { downloadImpl, isSmilesLookupTablePresent, saveSmilesLookupTable } from '../../Utility/Utils';
import { LoadingIndicatorDialog } from './LoadingIndicatorDialog';

export function selectPositions(dataset: Dataset, projection: IProjection) {
  const xChannel = projection.xChannel ?? 'x';
  const yChannel = projection.xChannel ?? 'y';
  return dataset.vectors.map((vector) => ({
    x: xChannel ? vector[xChannel] : 0,
    y: yChannel ? vector[yChannel] : 0,
  }));
}

function exportState(state: AppState) {
  const dump = UtilityActions.partialDump(state, [
    'dataset',
    'projectionColumns',
    'projects',
    'interfaceState',
    'lineUpInput',
    'pointColorMapping',
    'detailView',
  ]);
  // TODO: add lineup dump
  // before serialization, save the lineup dump in the pse dump
  // state.lineUpInput.lineup!.data.exportTable(lineUpInput.lineup!.data.getRankings()[0], {separator: ","}).then(x => downloadImpl(x, `lineup-export.csv`, 'application/csv'));
  // dump.lineupProviderDump = ...
  const dumpString = JSON.stringify(dump);
  downloadImpl(dumpString, `session.pse`, 'application/json');
}
function importState(files, dispatch) {
  if (files == null || files.length <= 0) {
    return;
  }
  const file = files[0];
  const fileName = file.name as string;

  if (fileName.endsWith('pse')) {
    const fileReader = new FileReader();
    fileReader.onload = (e) => {
      const dumpString = e.target.result.toString();
      dispatch(RootActions.hydrate(JSON.parse(dumpString)));
    };
    fileReader.readAsBinaryString(file);
  }
}

const mapStateToProps = (state: AppState) => ({});

const mapDispatchToProps = (dispatch) => ({
  setTriggerUpdate: (value) => dispatch(setTriggerUpdateAction(value)),
  resetViews: () => dispatch(CIME4RViewActions.resetViews()),
  hydrateState: (dump) => dispatch(RootActions.hydrate(dump)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  onDataSelected;
};

export const DatasetTabPanel = connector(({ onDataSelected, resetViews, setTriggerUpdate, hydrateState }: Props) => {
  const { clientConfig } = useVisynAppContext();
  const { cancellablePromise, cancelPromises } = useCancellablePromise();
  const abortController = new AbortController();
  const [refreshUploadedFiles, setRefreshUploadedFiles] = useState(0);
  const lookupFileInput = React.useRef<HTMLInputElement>();
  const [lookupUploadNote, setLookupUploadNote] = useState<string>(isSmilesLookupTablePresent());
  const stateFileInput = React.useRef<HTMLInputElement>();
  // const currentState = usePSESelector<AppState>((state) => state);
  const context = React.useContext(ReactReduxContext);
  const dispatch = useDispatch();

  const intermediateOnDataSelected = (dataset, state_dump?) => {
    resetViews();
    onDataSelected(dataset);
    if (state_dump != null) {
      hydrateState(state_dump);
    }
  };

  const triggerUpdate = (entry, state?: AppState) => {
    let stateDump: Partial<AppState> = null;
    if (state != null) {
      stateDump = UtilityActions.partialDump(state, [
        'dataset',
        'activeLine',
        'currentAggregation',
        'genericFingerprintAttributes',
        'handleDataset',
        'highlightedSequence',
        'hoverState',
        'mouseInteractionHooks',
        'projectionColumns',
        'projects',
        'interfaceState',
        'selectedLineBy',
        'viewTransform',
      ]);
    }
    new BackendCSVLoader().resolvePath(
      entry,
      (dataset) => {
        if (stateDump != null) {
          stateDump.multiples.multiples.ids.forEach((id) => {
            const multiple = stateDump.multiples.multiples.entities[id];
            const { workspace } = multiple.attributes;

            const projection = isEntityId(workspace) ? stateDump.multiples.projections.entities[workspace] : workspace;

            // Projection was not specified through x/y channels
            if (projection.positions) {
              projection.positions = selectPositions(dataset, projection);
            }
          });
        }

        dispatch(RootActions.loadDataset(dataset, stateDump));
        // intermediateOnDataSelected(dataset, stateDump);
      },
      cancellablePromise,
      null,
      abortController,
    );
  };

  React.useEffect(() => {
    setTriggerUpdate(triggerUpdate);
    // eslint-disable-next-line
  }, []);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      <UploadedFiles
        onChange={(entry) => {
          triggerUpdate(entry);
        }}
        refresh={refreshUploadedFiles}
      />

      <LoadingIndicatorDialog
        handleClose={() => {
          cancelPromises();
        }}
        area="global_loading_indicator"
      />

      {!clientConfig.publicVersion ? (
        <DatasetDrop
          onDatasetChange={(dataset) => {
            intermediateOnDataSelected(dataset);
            setRefreshUploadedFiles(refreshUploadedFiles + 1);
          }}
          cancellablePromise={cancellablePromise}
          abort_controller={abortController}
        />
      ) : null}

      {!clientConfig.publicVersion ? (
        <Box paddingLeft={2} paddingTop={2} paddingRight={2}>
          <Accordion variant="outlined" color="primary">
            <AccordionSummary expandIcon={<ExpandMore />}>
              <Typography variant="subtitle2">Advanced settings</Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Box>
                <Typography color="textSecondary" variant="body2">
                  Define a lookup table for molecule shortnames.{' '}
                  <Tooltip
                    title={
                      <Typography variant="subtitle2">
                        Upload a csv-file with the columns &#34;smiles&#34; and &#34;shortname&#34;. This mapping is then used to show a human readable name for
                        a molecule instead of the SMILES string.
                      </Typography>
                    }
                  >
                    <InfoOutlined fontSize="inherit" />
                  </Tooltip>
                </Typography>
                <input
                  style={{ display: 'none' }}
                  accept=".csv"
                  ref={lookupFileInput}
                  type="file"
                  onChange={(e) => {
                    saveSmilesLookupTable(e.target.files, (note) => {
                      setLookupUploadNote(() => note);
                    });
                  }}
                />
                <Button
                  fullWidth
                  variant="outlined"
                  aria-label="Define lookup table for shortnames of SMILES"
                  color="primary"
                  onClick={() => lookupFileInput?.current.click()}
                >
                  <ManageSearchIcon />
                  &nbsp;Select table
                </Button>
                <Typography variant="subtitle2" paddingTop={1} color="textSecondary">
                  {lookupUploadNote}
                </Typography>
              </Box>
              <Box paddingTop={2} paddingX={2}>
                <Divider orientation="horizontal" />
              </Box>
              <Box paddingTop={2}>
                <Typography variant="body2" paddingTop={1} color="textSecondary">
                  Save the current or load a previous session of the application.
                  <Tooltip
                    title={
                      <Typography variant="subtitle2">
                        You can download the current session of the application as a file and upload it again to restore the session. This can be useful if you
                        want to continue working at a later point in time or if you want to share your session.
                      </Typography>
                    }
                  >
                    <InfoOutlined fontSize="inherit" />
                  </Tooltip>
                </Typography>
                <Grid container>
                  <Grid item xs={6} paddingRight={1}>
                    <Button
                      fullWidth
                      variant="outlined"
                      onClick={() => {
                        exportState(context.store.getState());
                      }}
                    >
                      <GetApp />
                      &nbsp;Export
                    </Button>
                  </Grid>

                  <Grid item xs={6} paddingLeft={1}>
                    <input
                      style={{ display: 'none' }}
                      accept=".pse"
                      ref={stateFileInput}
                      // multiple
                      type="file"
                      onChange={(e) => {
                        importState(e.target.files, context.store.dispatch);
                      }}
                    />
                    <Button
                      fullWidth
                      variant="outlined"
                      onClick={() => {
                        stateFileInput?.current.click();
                      }}
                    >
                      <FileUpload />
                      &nbsp;Import
                    </Button>
                  </Grid>
                </Grid>
              </Box>
            </AccordionDetails>
          </Accordion>
        </Box>
      ) : null}
    </div>
  );
});
