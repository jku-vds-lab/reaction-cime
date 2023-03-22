import { Box, Divider, Grid, IconButton, Tooltip, Typography } from '@mui/material';
import { Dataset } from 'projection-space-explorer';
import React from 'react';
import DeleteIcon from '@mui/icons-material/Delete';
import { InfoOutlined } from '@mui/icons-material';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { formatLabel } from '../../Utility/Utils';
import { AppState } from '../../State/Store';

type Props = {
  dataset: Dataset;
  triggerDatasetUpdate;
  state: AppState;
};

export function ExceptionSettings({ dataset, triggerDatasetUpdate, state }: Props) {
  const [exceptions, setExceptions] = React.useState([]);

  const dropException = (index) => {
    const newEx = exceptions.filter((e, i) => i !== index);
    // setExceptions(new_ex)
    ReactionCIMEBackendFromEnv.updatePOIExceptions(dataset.info.path, newEx).then((res) => {
      if (res.msg !== 'ok') {
        alert(res.msg);
      }

      if (triggerDatasetUpdate != null) {
        triggerDatasetUpdate(
          {
            display: dataset.info.path,
            path: dataset.info.path,
            type: dataset.info.type,
            uploaded: true,
          },
          state,
        );
      }
    });
  };

  React.useEffect(() => {
    ReactionCIMEBackendFromEnv.loadPOIExceptions(dataset.info.path).then((res_constraints) => {
      setExceptions(res_constraints);
    });
  }, [dataset]);

  return (
    exceptions.length > 0 && (
      <div>
        <Box paddingY={2}>
          <Divider orientation="horizontal" />
        </Box>
        <Box paddingX={2} paddingTop={2} paddingBottom={1}>
          <Typography variant="subtitle2" gutterBottom>
            Exception settings
          </Typography>
          <Typography variant="body2" color="textSecondary" gutterBottom>
            Exceptions show all {state?.globalLabels.itemLabelPlural} close to a user-defined point in the scatterplot.{' '}
            <Tooltip
              title={
                <Typography variant="subtitle2">
                  To manually define exceptions, right-click on any coordinates in the scatter to open the context menu and select &quot;Show{' '}
                  {state?.globalLabels.itemLabelPlural} around this point&quot;. Exceptions are stored persistently and will therefore be available in future
                  sessions.
                </Typography>
              }
            >
              <InfoOutlined fontSize="inherit" />
            </Tooltip>
          </Typography>
        </Box>
        {exceptions.map((exc, i) => (
          // eslint-disable-next-line react/no-array-index-key
          <Grid key={`exception${i}`} container paddingTop={0}>
            <Grid item xs={3} textAlign="right">
              <Tooltip title={<Typography variant="subtitle2">Remove this exception.</Typography>}>
                <IconButton onClick={() => dropException(i)}>
                  <DeleteIcon fontSize="large" />
                </IconButton>
              </Tooltip>
            </Grid>
            <Grid item xs={9} margin="auto">
              <Tooltip placement="right" title={<Typography variant="subtitle2">Coordinates of the exception point.</Typography>}>
                <Typography width="fit-content">
                  {exc.x_col}: {formatLabel(exc.x_coord)} | {exc.y_col}: {formatLabel(exc.y_coord)}
                </Typography>
              </Tooltip>
            </Grid>
          </Grid>
        ))}
      </div>
    )
  );
}
