import * as React from 'react';
import { Box, Button, Tooltip, Typography } from '@mui/material';
import { Dataset, IProjection, RootActions, useCancellablePromise, UtilityActions } from 'projection-space-explorer';
import { connect, ConnectedProps } from 'react-redux';
import { useState } from 'react';
import ManageSearchIcon from '@mui/icons-material/ManageSearch';
import { DatasetDrop } from './DatasetDrop';
import { AppState, CIME4RViewActions } from '../../State/Store';
import { UploadedFiles } from './UploadedFiles';
import { BackendCSVLoader } from './BackendCSVLoader';
import { setTriggerUpdate as setTriggerUpdateAction } from '../../State/HandleDatasetDuck';
import { saveSmilesLookupTable } from '../../Utility/Utils';
import { LoadingIndicatorDialog } from './LoadingIndicatorDialog';

function selectPositions(dataset: Dataset, projection: IProjection) {
  const xChannel = projection.xChannel ?? 'x';
  const yChannel = projection.xChannel ?? 'y';
  return dataset.vectors.map((vector) => ({
    x: xChannel ? vector[xChannel] : 0,
    y: yChannel ? vector[yChannel] : 0,
  }));
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
  const { cancellablePromise, cancelPromises } = useCancellablePromise();
  const abortController = new AbortController();
  const [refreshUploadedFiles, setRefreshUploadedFiles] = useState(0);
  const lookupFileInput = React.useRef<any>();

  const intermediateOnDataSelected = (dataset, state_dump?) => {
    resetViews();
    onDataSelected(dataset);
    if (state_dump != null) {
      hydrateState(state_dump);
    }
  };

  const triggerUpdate = (entry, state?: AppState) => {
    let stateDump = null;
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
      ]);
    }
    new BackendCSVLoader().resolvePath(
      entry,
      (dataset) => {
        if (stateDump != null) {
          // we have to update the workspace positions manually to the new positions
          const newProjectionEntities = { ...stateDump.multiples.projections.entities };
          stateDump.multiples.multiples.ids.forEach((id) => {
            const active = stateDump.multiples.multiples.entities[id];
            const workspaceId = active.attributes.workspace;
            const workspace = state.multiples.multiples.entities[id].attributes.workspace as IProjection;
            const newPositions = selectPositions(dataset, workspace);
            const newWorkspaceIdPosition = { ...newProjectionEntities[workspaceId] };
            newWorkspaceIdPosition.positions = newPositions;
            newProjectionEntities[workspaceId] = newWorkspaceIdPosition;
          });
          const newProjections = { ...stateDump.multiples.projections, entities: newProjectionEntities };
          const newMultiples = { ...stateDump.multiples, projections: newProjections };
          stateDump = { ...stateDump, multiples: newMultiples };
        }
        intermediateOnDataSelected(dataset, stateDump);
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
      <Box paddingLeft={2} paddingTop={2}>
        <Typography variant="subtitle2" gutterBottom>
          Custom datasets (drag and drop)
        </Typography>
      </Box>

      <DatasetDrop
        onDatasetChange={(dataset) => {
          intermediateOnDataSelected(dataset);
          setRefreshUploadedFiles(refreshUploadedFiles + 1);
        }}
        cancellablePromise={cancellablePromise}
        abort_controller={abortController}
      />

      <Box paddingLeft={2} paddingTop={2}>
        <Typography variant="subtitle2" gutterBottom>
          Predefined datasets
        </Typography>
      </Box>

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
      <input
        style={{ display: 'none' }}
        accept=".csv"
        ref={lookupFileInput}
        type="file"
        onChange={(e) => {
          saveSmilesLookupTable(e.target.files);
        }}
      />
      <Box paddingLeft={2} paddingTop={2} paddingRight={2}>
        <Tooltip title={'Select a lookup table for shortnames of molecules. It has to be a csv-file with the columns "smiles" and "shortname".'}>
          <Button
            fullWidth
            variant="outlined"
            aria-label="Define lookup table for shortnames of SMILES"
            color="primary"
            onClick={() => lookupFileInput.current.click()}
          >
            <ManageSearchIcon />
            &nbsp;Define lookup table
          </Button>
        </Tooltip>
      </Box>
    </div>
  );
});
