import { Box, Divider, Grid, IconButton, Tooltip, Typography } from '@mui/material';
import { Dataset } from 'projection-space-explorer';
import React from 'react';
import DeleteIcon from '@mui/icons-material/Delete';
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { formatLabel } from '../../Utility/Utils';

type Props = {
  dataset: Dataset;
  triggerDatasetUpdate;
  state;
};

export function ExceptionSettings({dataset, triggerDatasetUpdate, state}:Props) {
    const [exceptions, setExceptions] = React.useState([]);

    const dropException = (index, state) => {
        const new_ex = exceptions.filter((e, i) => i !== index)
        // setExceptions(new_ex)
        ReactionCIMEBackendFromEnv.updatePOIExceptions(dataset.info.path, new_ex).then((res) => {
            if(res.msg !== "ok")
                alert(res.msg)

            if(triggerDatasetUpdate != null){
                triggerDatasetUpdate({
                    display: dataset.info.path,
                    path: dataset.info.path,
                    type: dataset.info.type,
                    uploaded: true
                }, state)
            }
        })
    }

    React.useEffect(() => {
        ReactionCIMEBackendFromEnv.loadPOIExceptions(dataset.info.path).then((res_constraints) => {
            setExceptions(res_constraints)
        })
    }, [dataset])

    return exceptions.length > 0 && <div>
            <Box paddingY={2}>
                <Divider orientation="horizontal" />
            </Box>
            <Box paddingLeft={2} paddingRight={2}>
                <Typography variant="subtitle2" gutterBottom>
                    Define Exceptions
                </Typography>
            </Box>
            {exceptions.map((exc, i) => 
            <Grid key={"exception"+i} container paddingTop={0}>
                <Grid item xs={3} textAlign={"right"}>
                    <Tooltip title={"Remove exception"}>
                        <IconButton onClick={() => dropException(i, state)}><DeleteIcon fontSize="large"/></IconButton>
                    </Tooltip>
                </Grid>
                <Grid item xs={9} margin={"auto"}>
                    <Typography>
                        {exc.x_col}: {formatLabel(exc.x_coord)} | {exc.y_col}: {formatLabel(exc.y_coord)}
                    </Typography>
                </Grid>
            </Grid>
            )}
        </div>
}
};
