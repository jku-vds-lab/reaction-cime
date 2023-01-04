import * as React from 'react';
import { Grid } from '@mui/material';
import { DragAndDrop } from 'projection-space-explorer';
import { BackendCSVLoader } from './BackendCSVLoader';

export function DatasetDrop({ onDatasetChange, cancellablePromise, abort_controller }) {
  return (
    <Grid container item alignItems="stretch" justifyContent="center" direction="column" style={{ padding: '16px' }}>
      {/** @ts-ignore* */}
      <DragAndDrop
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
      </DragAndDrop>
    </Grid>
  );
}
