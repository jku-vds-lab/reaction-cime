import {
  Box,
  Button,
  CircularProgress,
  CircularProgressProps,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  IconButton,
  List,
  ListItem,
  ListItemButton,
  ListItemText,
  Tooltip,
  Typography,
} from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';
import { useCancellablePromise, usePSESelector } from 'projection-space-explorer';
import React, { CSSProperties } from 'react';
import RefreshIcon from '@mui/icons-material/Refresh';
import { useVisynAppContext } from 'visyn_core/app';
import { userSession } from 'visyn_core/security';
import { useDispatch, useSelector } from 'react-redux';

import { AppState } from '../../State/Store';
import { deleteProject, syncProjects } from '../../State/ProjectsDuck';

const textOverflowStyle: CSSProperties = {
  textOverflow: 'ellipsis',
  overflow: 'hidden',
  whiteSpace: 'nowrap',
};

function CircularProgressWithLabel(props: CircularProgressProps & { value: number }) {
  return (
    <Box sx={{ position: 'relative', display: 'inline-flex' }}>
      <CircularProgress />
      <Box
        sx={{
          top: 0,
          left: 0,
          bottom: 0,
          right: 0,
          position: 'absolute',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
        }}
      >
        <Typography display="none" variant="caption" component="div" color="text.secondary">{`${Math.round(props.value)}`}</Typography>
      </Box>
    </Box>
  );
}

export function UploadedFiles({ onChange, refresh }) {
  const { clientConfig } = useVisynAppContext();
  const dispatch = useDispatch();
  const files = useSelector((state: AppState) => Object.values(state.projects.projects.entities));
  const { cancellablePromise } = useCancellablePromise();
  const dataset = usePSESelector((state) => state.dataset);
  const [deleteDialog, setDeleteDialog] = React.useState(false);
  const [item, setItem] = React.useState<{ name: string; id: string }>();

  const updateFiles = () => {
    dispatch(syncProjects());
  };

  const updateRef = React.useRef(updateFiles);
  updateRef.current = updateFiles;

  React.useEffect(() => {
    if (files?.some((file) => file.file_status.startsWith('Processing'))) {
      const interval = setInterval(() => {
        updateRef.current();
      }, 5000);

      return () => clearInterval(interval);
    }
    return undefined;
  }, [files]);

  React.useEffect(() => {
    updateFiles();
    // eslint-disable-next-line
  }, [refresh]);

  return (
    files && (
      <div>
        <Box paddingLeft={2} paddingRight={2} paddingTop={2}>
          {!clientConfig.publicVersion && (
            <Typography variant="subtitle2" gutterBottom>
              Select dataset{' '}
              <Tooltip title="Refresh dataset list">
                <Button onClick={() => updateFiles()}>
                  <RefreshIcon style={{ fontSize: '1.25rem' }} />
                </Button>
              </Tooltip>
            </Typography>
          )}
          <List
            subheader={<li />}
            style={{ backgroundColor: 'white', border: '1px solid lightgrey', borderRadius: '4px', overflowY: 'auto', maxHeight: '400px' }}
          >
            {files.length === 0 ? <Typography color="gray">No datasets uploaded</Typography> : null}
            {files.map((file) => (
              <ListItem
                disablePadding
                key={file.id}
                selected={file.id === dataset?.info?.path}
                secondaryAction={
                  file.file_status.startsWith('Processing') ? (
                    <CircularProgressWithLabel value={Number.parseInt(file.file_status.substring('Processing_'.length), 10)} />
                  ) : !clientConfig.publicVersion && userSession.canWrite(file) ? (
                    <Tooltip placement="right" title={<Typography variant="subtitle2">Permanently delete dataset &quot;{file.name}&quot;.</Typography>}>
                      <IconButton
                        edge="end"
                        aria-label="delete"
                        onClick={(event) => {
                          event.preventDefault();
                          setDeleteDialog(true);
                          setItem(file);
                        }}
                      >
                        <DeleteIcon />
                      </IconButton>
                    </Tooltip>
                  ) : (
                    <div />
                  )
                }
              >
                <Tooltip placement="right" title={<Typography variant="subtitle2">Load dataset &quot;{file.name}&quot;.</Typography>}>
                  <ListItemButton
                    style={{ width: '100%' }}
                    key={file.id}
                    disabled={file.file_status.startsWith('Error') || file.file_status.startsWith('Processing')}
                    data-cy="uploaded-data-list-item"
                    href={`/?project=${file.id}`}
                    component="a"
                    target="_self"
                  >
                    <ListItemText
                      primary={file.name}
                      secondary={file.file_status.startsWith('Error') ? file.file_status : `By ${file.creator}`}
                      primaryTypographyProps={{
                        style: textOverflowStyle,
                      }}
                      secondaryTypographyProps={{
                        style: textOverflowStyle,
                        color: file.file_status.startsWith('Error') ? 'error' : 'textSecondary',
                      }}
                    />
                  </ListItemButton>
                </Tooltip>
              </ListItem>
            ))}
          </List>
        </Box>
        <Dialog open={deleteDialog} onClose={() => setDeleteDialog(false)} aria-labelledby="alert-dialog-title" aria-describedby="alert-dialog-description">
          <DialogTitle id="alert-dialog-title">Delete dataset &quot;{item?.name}&quot;</DialogTitle>
          <DialogContent>
            <DialogContentText id="alert-dialog-description">
              Are you sure you want to delete dataset &quot;{item?.name}&quot;? This action cannot be undone.
            </DialogContentText>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setDeleteDialog(false)}>Cancel</Button>
            <Button
              onClick={() => {
                dispatch(deleteProject(item?.id));
                setDeleteDialog(false);
              }}
              autoFocus
            >
              Delete
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  );
}
