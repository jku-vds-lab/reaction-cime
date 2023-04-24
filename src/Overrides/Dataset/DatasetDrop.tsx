import * as React from 'react';
import { Alert, Button, Grid, IconButton, Snackbar } from '@mui/material';
import { Dataset } from 'projection-space-explorer';
import CloseIcon from '@mui/icons-material/Close';
import { useDispatch } from 'react-redux';
import { BackendCSVLoader } from './BackendCSVLoader';
import { Project, ProjectActions } from '../../State/ProjectsDuck';

export function DatasetDrop({ onDatasetChange, cancellablePromise, abort_controller }) {
  const fileInput = React.useRef<HTMLInputElement>();
  const [msg, setMsg] = React.useState('');

  const dispatch = useDispatch();

  const openSnack = (value: string) => {
    setMsg(value);
  };

  const handleClose = (event: React.SyntheticEvent | Event, reason?: string) => {
    if (reason === 'clickaway') {
      return;
    }

    setMsg('');
  };

  const action = (
    <IconButton size="small" aria-label="close" color="inherit" onClick={handleClose}>
      <CloseIcon fontSize="small" />
    </IconButton>
  );

  return (
    <Grid container item alignItems="stretch" justifyContent="center" direction="column" style={{ padding: '16px' }}>
      <input
        style={{ display: 'none' }}
        accept=".csv,.zip"
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

          if (fileName.endsWith('csv') || fileName.endsWith('zip')) {
            const loader = new BackendCSVLoader();
            loader.upload_file(
              file,
              (project: Project) => {
                dispatch(ProjectActions.addProject(project));
              },
              cancellablePromise,
              abort_controller,
            );
            e.target.value = null;
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
      <Snackbar open={msg !== ''} autoHideDuration={10000} onClose={handleClose} message={msg} action={action}>
        <Alert severity="info">{msg}</Alert>
      </Snackbar>
    </Grid>
  );
}
