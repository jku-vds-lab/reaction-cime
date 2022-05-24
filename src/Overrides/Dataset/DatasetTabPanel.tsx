import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    Tooltip,
    Typography,
  } from "@mui/material";
import { RootActions, useCancellablePromise, UtilityActions } from "projection-space-explorer";
import { usePromiseTracker } from "react-promise-tracker";
import { DatasetDrop } from "./DatasetDrop";
import Loader from "react-loader-spinner";
import { AppState, CIME4RViewActions } from "../../State/Store";
import { connect, ConnectedProps } from "react-redux";
import { UploadedFiles } from "./UploadedFiles";
import React, { useState } from "react";
import { BackendCSVLoader } from "./BackendCSVLoader";
import { setTriggerUpdate } from "../../State/HandleDatasetDuck";
import { save_smiles_lookup_table } from "../../Utility/Utils";
import ManageSearchIcon from '@mui/icons-material/ManageSearch';

  export const LoadingIndicatorView = (props) => {
    const { promiseInProgress } = usePromiseTracker({ area: props.area });
    
    return (
      promiseInProgress && (
        <div
          style={{
            width: "100%",
            height: "100",
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
          }}
        >
          <Loader type="ThreeDots" color="#2BAD60" height="100" width="100" />
        </div>
      )
    );
  };
  
  export const LoadingIndicatorDialog = (props) => {
    const { promiseInProgress } = usePromiseTracker({ area: props.area });
  
    return (
      <Dialog maxWidth="lg" open={promiseInProgress}>
        {" "}
        {/*onClose={props.handleClose} */}
        <DialogContent>
          <LoadingIndicatorView area={props.area} />
        </DialogContent>
        <DialogActions>
          <Button onClick={props.handleClose}>Cancel</Button>
        </DialogActions>
      </Dialog>
    );
  };



  const mapStateToProps = (state: AppState) => ({
  });
  
  const mapDispatchToProps = (dispatch) => ({
    setTriggerUpdate: value => dispatch(setTriggerUpdate(value)),
    resetViews: () => dispatch(CIME4RViewActions.resetViews()),
    hydrateState: (dump) => dispatch(RootActions.hydrate(dump)),
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
    onDataSelected
  };

  export const DatasetTabPanel = connector(({onDataSelected, resetViews, setTriggerUpdate, hydrateState}: Props) => {
    const { cancellablePromise, cancelPromises } = useCancellablePromise();
    let abort_controller = new AbortController();
    const [refreshUploadedFiles, setRefreshUploadedFiles] = useState(0);
    let lookupFileInput = React.useRef<any>();

    const intermediateOnDataSelected = (dataset, state_dump?) => {
      resetViews();
      onDataSelected(dataset);
      if(state_dump != null){
        console.log("hydrate state")
        hydrateState(state_dump)
      }
      
    }

    const triggerUpdate = (entry, state?) => {
      let state_dump = null;
      if(state != null)
        state_dump = UtilityActions.partialDump(state, ['dataset', 'activeLine', 'currentAggregation', 'genericFingerprintAttributes', 'handleDataset', 'highlightedSequence', 'hoverState', 'mouseInteractionHooks', 'projectionColumns', 'projects', 'interfaceState', 'selectedLineBy']);

      new BackendCSVLoader().resolvePath(
        entry,
        (dataset) => {
          intermediateOnDataSelected(dataset, state_dump);
        },
        cancellablePromise,
        null,
        abort_controller
      );
    }

    React.useEffect(() => {
      setTriggerUpdate(triggerUpdate)
      // eslint-disable-next-line
    }, [])
  
    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={2}>
          <Typography variant="subtitle2" gutterBottom>
            {"Custom Datasets (Drag and Drop)"}
          </Typography>
        </Box>
  
        <DatasetDrop
          onDatasetChange={(dataset) => {
            intermediateOnDataSelected(dataset);
            setRefreshUploadedFiles(refreshUploadedFiles + 1);
          }}
          cancellablePromise={cancellablePromise}
          abort_controller={abort_controller}
        />
  
      <Box paddingLeft={2} paddingTop={2}>
        <Typography variant="subtitle2" gutterBottom>
          {"Predefined Datasets"}
        </Typography>
      </Box>

      <UploadedFiles
        onChange={(entry) => {
          triggerUpdate(entry)
        }}
        refresh={refreshUploadedFiles}
      />

      <LoadingIndicatorDialog
        handleClose={() => {
          cancelPromises();
        }}
        area={"global_loading_indicator"}
      />
      <input
          style={{ display: 'none' }}
          accept={".csv"}
          ref={lookupFileInput}
          type="file"
          onChange={(e) => {
            save_smiles_lookup_table(e.target.files);
          }}
      />
      <Box paddingLeft={2} paddingTop={2} paddingRight={2}>
        <Tooltip title={'Select a lookup table for shortnames of molecules. It has to be a csv-file with the columns "smiles" and "shortname".'}>
          <Button fullWidth variant="outlined" aria-label="Define lookup table for shortnames of SMILES" color="primary" onClick={() => lookupFileInput.current.click()}>
            <ManageSearchIcon />&nbsp;Define lookup table
          </Button>
        </Tooltip>
      </Box>

      </div>
    );
  });