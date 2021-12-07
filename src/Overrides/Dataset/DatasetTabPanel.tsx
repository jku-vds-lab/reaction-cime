import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    Typography,
  } from "@mui/material";
import { useCancellablePromise } from "projection-space-explorer";
import { usePromiseTracker } from "react-promise-tracker";
import { DatasetDrop } from "./DatasetDrop";
import Loader from "react-loader-spinner";
import { AppState } from "../../State/Store";
import { connect, ConnectedProps } from "react-redux";
import { setAggregateDatasetAction } from "../../State/AggregateDatasetDuck";
import { AggregateDataset } from "../AggregationTabPanel/AggregateDataset";
import { UploadedFiles } from "./UploadedFiles";
import { useState } from "react";
import { BackendCSVLoader } from "./BackendCSVLoader";
  
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
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
    onDataSelected
  };

  export const DatasetTabPanel = connector(({onDataSelected}: Props) => {
    const { cancellablePromise, cancelPromises } = useCancellablePromise();
    let abort_controller = new AbortController();
    const [refreshUploadedFiles, setRefreshUploadedFiles] = useState(0);
  
    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={2}>
          <Typography variant="subtitle2" gutterBottom>
            {"Custom Datasets (Drag and Drop)"}
          </Typography>
        </Box>
  
        <DatasetDrop
          onDatasetChange={(dataset) => {
            onDataSelected(dataset);
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
          // setEntry(entry);
          new BackendCSVLoader().resolvePath(
            entry,
            (dataset) => {
              onDataSelected(dataset);
            },
            cancellablePromise,
            null,
            abort_controller
        );
        }}
        refresh={refreshUploadedFiles}
      />

      <LoadingIndicatorDialog
        handleClose={() => {
          cancelPromises();
        }}
        area={"global_loading_indicator"}
      />
      </div>
    );
  });