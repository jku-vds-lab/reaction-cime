import { Box, Button, IconButton, List, ListItem, ListItemButton, ListItemText, Tooltip, Typography } from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';
import { useCancellablePromise } from 'projection-space-explorer';
import React, { CSSProperties } from 'react';
import RefreshIcon from '@mui/icons-material/Refresh';
import { trackPromise } from 'react-promise-tracker';
import { ISecureItem, userSession, useVisynAppContext } from 'visyn_core';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { LoadingIndicatorView } from './LoadingIndicatorDialog';

const textOverflowStyle: CSSProperties = {
  textOverflow: 'ellipsis',
  overflow: 'hidden',
  whiteSpace: 'nowrap',
};

const loadingArea = 'update_uploaded_files_list';
export function UploadedFiles({ onChange, refresh }) {
  const { clientConfig } = useVisynAppContext();
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
        <Box paddingLeft={2} paddingRight={2} paddingTop={2}>
          {clientConfig.publicVersion && (
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
            <LoadingIndicatorView area={loadingArea} />

            {files.map((file) => (
              <ListItem
                disablePadding
                key={file.id}
                secondaryAction={
                  !clientConfig.publicVersion && userSession.canWrite(file) ? (
                    <Tooltip placement="right" title={<Typography variant="subtitle2">Permanently delete dataset &quot;{file.name}&quot;</Typography>}>
                      <IconButton
                        edge="end"
                        aria-label="delete"
                        onClick={(event) => {
                          event.preventDefault();
                          handleDelete(file.id);
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
                <Tooltip placement="right" title={<Typography variant="subtitle2">Load dataset &quot;{file.name}&quot;</Typography>}>
                  <ListItemButton
                    style={{ width: '100%' }}
                    key={file.id}
                    data-cy="uploaded-data-list-item"
                    href={`/?project=${file.id}`}
                    component="a"
                    target="_self"
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
                  </ListItemButton>
                </Tooltip>
              </ListItem>
            ))}
          </List>
        </Box>
      </div>
    )
  );
}
