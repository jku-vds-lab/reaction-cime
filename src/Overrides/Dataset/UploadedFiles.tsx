import { Button, Grid, IconButton, List, ListItem, ListItemSecondaryAction, ListItemText, ListSubheader } from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';
import { DatasetType, useCancellablePromise } from 'projection-space-explorer';
import React, { CSSProperties } from 'react';
import RefreshIcon from '@mui/icons-material/Refresh';
import { trackPromise } from 'react-promise-tracker';
import { ISecureItem, userSession } from 'visyn_core';
import { DEMO } from '../../constants';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { LoadingIndicatorView } from './LoadingIndicatorDialog';

const textOverflowStyle = {
  textOverflow: 'ellipsis',
  overflow: 'hidden',
  whiteSpace: 'nowrap',
} satisfies CSSProperties;

const loadingArea = 'update_uploaded_files_list';
export function UploadedFiles({ onChange, refresh }) {
  const [files, setFiles] = React.useState<({ name: string; id: string } & ISecureItem)[]>([]);
  const { cancellablePromise } = useCancellablePromise();

  const updateFiles = () => {
    trackPromise(
      cancellablePromise(ReactionCIMEBackendFromEnv.getUploadedFiles())
        .then((data) => {
          if (data) {
            setFiles(data);
          }
        })
        .catch((error) => console.log(error)),
      loadingArea,
    );
  };

  React.useEffect(() => {
    updateFiles();
    // eslint-disable-next-line
  }, [refresh]);

  const handleClick = (entry) => {
    onChange(entry);
  };

  const handleDelete = (file: string) => {
    cancellablePromise(ReactionCIMEBackendFromEnv.deleteFile(file))
      .then((response) => {
        if (response && response.deleted === 'true') setFiles(files.filter((f) => f.id !== file));
      })
      .catch((error) => console.log(error));
  };

  return (
    files && (
      <div>
        <Grid item style={{ overflowY: 'auto', flex: '1 1 auto', maxHeight: '400px' }}>
          <List subheader={<li />} style={{ backgroundColor: 'white' }}>
            {!DEMO && (
              <ListSubheader>
                Uploaded files{' '}
                <Button onClick={() => updateFiles()}>
                  <RefreshIcon style={{ fontSize: '1.25rem' }} />
                </Button>
              </ListSubheader>
            )}
            {DEMO && <ListSubheader>Select Dataset</ListSubheader>}
            {files.map((file) => (
              <ListItem
                key={file.id}
                data-cy="uploaded-data-list-item"
                button
                onClick={() => {
                  handleClick({
                    display: file.name,
                    path: file.id,
                    type: DatasetType.Chem,
                    uploaded: true, // indicates that file is already uploaded
                  });
                }}
              >
                <ListItemText
                  primary={file.name}
                  secondary={`By ${file.creator}`}
                  primaryTypographyProps={{
                    style: textOverflowStyle,
                  }}
                  secondaryTypographyProps={{
                    style: textOverflowStyle,
                  }}
                />
                {!DEMO && userSession.canWrite(file) && (
                  <ListItemSecondaryAction
                    onClick={() => {
                      handleDelete(file.id);
                    }}
                  >
                    <IconButton edge="end" aria-label="delete">
                      <DeleteIcon />
                    </IconButton>
                  </ListItemSecondaryAction>
                )}
              </ListItem>
            ))}
          </List>
          <LoadingIndicatorView area={loadingArea} />
        </Grid>
      </div>
    )
  );
}
