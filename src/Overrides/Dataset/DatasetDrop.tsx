import { Grid } from "@mui/material";
import { DragAndDrop } from "projection-space-explorer";
import { BackendCSVLoader } from "./BackendCSVLoader";

export const DatasetDrop = ({
  onDatasetChange,
  cancellablePromise,
  abort_controller,
}) => {

  return (
    <Grid
      container
      item
      alignItems="stretch"
      justifyContent="center"
      direction="column"
      style={{ padding: "16px" }}
    >
      {/**@ts-ignore**/}
      <DragAndDrop
        accept=".csv"
        handleDrop={(files) => {
          if (files == null || files.length <= 0) {
            return;
          }

          var file = files[0];
          var fileName = file.name as string;

          if (fileName.endsWith("csv")) {
            abort_controller = new AbortController();
            new BackendCSVLoader().resolveContent(
                file,
                onDatasetChange,
                cancellablePromise,
                abort_controller
            );
          }
        }}
      >
        <div style={{ height: 200 }}></div>
      </DragAndDrop>
    </Grid>
  );
};
