import * as React from 'react';
import { TextField, Typography } from '@mui/material';
import { Box } from '@mui/system';

type Props = {
  title: string;
  target: string;
  setRange: any;
  range: { min: number; max: number };
};

export function MinMaxNumberInput({ title, target, setRange, range }: Props) {
  return (
    <Box paddingTop={1}>
      <Typography variant="subtitle2" gutterBottom>
        {title}
      </Typography>
      <TextField
        title={`range minimum for ${target}`}
        type="number"
        id="knn-textfield"
        label="min"
        variant="outlined"
        onChange={(event) => {
          const newVal = parseFloat(event.target.value);
          const current = { ...range };
          if (!Number.isNaN(newVal) && newVal < range.max) {
            current.min = newVal;
            setRange(current);
          } else {
            current.min = current.max;
            setRange(current);
          }
        }}
        style={{ width: '50%' }}
        size="small"
        value={range.min}
      />
      <TextField
        title={`range maximum for ${target}`}
        type="number"
        id="knn-textfield"
        label="max"
        variant="outlined"
        onChange={(event) => {
          const newVal = parseFloat(event.target.value);
          const current = { ...range };
          if (!Number.isNaN(newVal) && newVal > range.min) {
            current.max = newVal;
            setRange(current);
          } else {
            current.max = current.min;
            setRange(current);
          }
        }}
        style={{ width: '50%' }}
        size="small"
        value={range.max}
      />
    </Box>
  );
}
