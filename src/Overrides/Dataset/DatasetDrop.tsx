import * as React from 'react';
import { Button, Grid } from '@mui/material';
import { BackendCSVLoader } from './BackendCSVLoader';

export function DatasetDrop({ onDatasetChange, cancellablePromise, abort_controller }) {
  const fileInput = React.useRef<HTMLInputElement>();

  return (
    <Grid container item alignItems="stretch" justifyContent="center" direction="column" style={{ padding: '16px' }}>
      <input
        style={{ display: 'none' }}
        accept=".csv"
        ref={fileInput}
        data-cy="upload-file-input"
        // multiple
        type="file"
        onChange={(e) => {
          const { files } = e.target;
          if (files == null || files.length <= 0) {
            return;
          }

          const file = files[0];
          const fileName = file.name as string;

          if (fileName.endsWith('csv')) {
            // abort_controller = new AbortController();
            new BackendCSVLoader().resolveContent(file, onDatasetChange, cancellablePromise, abort_controller);
          }
        }}
      />
      <Button
        variant="outlined"
        component="span"
        onClick={() => {
          const fi = fileInput.current;
          fi.click();
        }}
      >
        Upload new dataset
      </Button>
      {/* <DragAndDrop
        accept=".csv"
        handleDrop={(files) => {
          if (files == null || files.length <= 0) {
            return;
          }

          const file = files[0];
          const fileName = file.name as string;

          if (fileName.endsWith('csv')) {
            // abort_controller = new AbortController();
            new BackendCSVLoader().resolveContent(file, onDatasetChange, cancellablePromise, abort_controller);
          }
        }}
      >
        <div style={{ height: 200 }} />
      </DragAndDrop> */}
    </Grid>
  );
}
