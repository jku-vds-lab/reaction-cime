import { Grid, IconButton, ToggleButton, Tooltip, Typography } from '@mui/material';
import React from 'react';
import DeleteIcon from '@mui/icons-material/Delete';
import { Dataset } from 'projection-space-explorer';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { mapSmilesToShortname } from '../../Utility/Utils';

type Props = {
  col: string;
  value: number[];
  setValue: (value: number[]) => void;
  remove: (col: string) => void;
  dataset: Dataset;
};

export function CategoryFilter({ col, value, setValue, remove, dataset }: Props) {
  const [catValues, setCatValues] = React.useState([]);

  React.useEffect(() => {
    ReactionCIMEBackendFromEnv.loadCategoryValues(dataset.info.path, col).then((result) => {
      setCatValues(result.values);
    });
  }, [col, dataset.info.path]);

  return (
    <Grid container paddingTop={0}>
      <Grid item xs={3} textAlign="right">
        <Tooltip title={<Typography variant="subtitle2">Remove {col} filter.</Typography>}>
          <IconButton
            onClick={() => {
              remove(col);
            }}
          >
            <DeleteIcon fontSize="large" />
          </IconButton>
        </Tooltip>
      </Grid>
      <Grid item xs={9}>
        <Typography id={`filter_${col}`} marginBottom="0px" style={{ textOverflow: 'ellipsis', overflow: 'hidden' }}>
          {col}
        </Typography>
        <div>
          {catValues.map((val) => {
            return (
              <Tooltip placement="right" title={<Typography variant="subtitle2">{val}</Typography>} key={val}>
                <ToggleButton
                  // fullWidth
                  style={{ margin: '1px', textOverflow: 'ellipsis', overflow: 'hidden' }}
                  color="primary"
                  value={val}
                  selected={value.includes(val)}
                  onChange={() => {
                    if (value.includes(val)) {
                      // if previously checked, remove from list
                      setValue([...new Set(value.filter((v) => v !== val))]);
                    } else {
                      // if previously unchecked, add to list
                      setValue([...new Set([...value, val])]);
                    }
                  }}
                >
                  &nbsp;{dataset.columns[col].metaInformation.imgSmiles ? mapSmilesToShortname(val) : val} (
                  {dataset.vectors.filter((vector) => vector[col] === val).length})&nbsp;
                </ToggleButton>
              </Tooltip>
            );
          })}
        </div>
      </Grid>
    </Grid>
  );
}
