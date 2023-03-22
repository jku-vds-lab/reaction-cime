/* eslint-disable guard-for-in */
import { Box, Button, Tooltip, Typography } from '@mui/material';
import { Dataset } from 'projection-space-explorer';
import React from 'react';
import FilterAltIcon from '@mui/icons-material/FilterAlt';
import RotateLeftIcon from '@mui/icons-material/RotateLeft';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { CategoryFilter } from './CategoryFilter';
import { RangeFilter } from './RangeFilter';
import { AppState } from '../../State/Store';

export const updateBackendConstraints = (dimensions: Record<string, any>, dataset, triggerDatasetUpdate, state) => {
  const constraintDimensions = dimensions;
  const allConstraints = [];
  for (const i in constraintDimensions) {
    const constDimension = constraintDimensions[i];
    if (constDimension.isNum) {
      const constraintObject = { col: i, operator: 'BETWEEN', val1: constDimension.val[0], val2: constDimension.val[1] };
      allConstraints.push(constraintObject);
    } else {
      const constraintarray = constDimension.val;
      for (const j in constraintarray) {
        const constraintObject = { col: i, operator: 'EQUALS', val1: constraintarray[j], val2: constraintarray[j] };
        allConstraints.push(constraintObject);
      }
    }
  }

  ReactionCIMEBackendFromEnv.updatePOIConstraints(dataset.info.path, allConstraints).then((res) => {
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

type Props = {
  dataset: Dataset;
  removeFilter: (col) => void;
  constraintCols: string[];
  constraints: { col: string; operator: string; val1: string; val2: string }[];
  triggerDatasetUpdate;
  state: AppState;
};

export function FilterSettings({ dataset, removeFilter, constraintCols, constraints, triggerDatasetUpdate, state }: Props) {
  const [filterValues, setFilterValues] = React.useState({});

  React.useEffect(() => {
    // initialize filters
    if (constraints != null) {
      const tempFilterValues = { ...filterValues };
      constraintCols.forEach((col) => {
        const currentConstraints = constraints.filter((item) => item.col === col);

        let between = [undefined, undefined];
        const equals = [];
        if (currentConstraints.length > 0) {
          currentConstraints.forEach((item) => {
            const con = item;
            if (con.operator === 'BETWEEN') {
              // if there are several between operators, we choose the minimal minimum and the maximal maximum, because we can only handle one range
              if (between[0] == null || between[1] == null) {
                between = [parseFloat(con.val1), parseFloat(con.val2)];
              } else {
                between[0] = Math.min(parseFloat(con.val1), between[0]);
                between[1] = Math.max(parseFloat(con.val2), between[1]);
              }
            } else if (con.operator === 'EQUALS') {
              equals.push(con.val1);
            }
          });
        }

        if (dataset.columns[col].isNumeric) {
          tempFilterValues[col] = { isNum: true, val: between };
        } else {
          tempFilterValues[col] = { isNum: false, val: equals };
        }
      });

      setFilterValues(tempFilterValues);
    }
    // eslint-disable-next-line
  }, [constraints]);

  React.useEffect(() => {
    if (constraintCols != null) {
      const tempFilterValues = { ...filterValues };
      constraintCols.forEach((col) => {
        if (!(col in tempFilterValues)) {
          if (dataset.columns[col].isNumeric) {
            tempFilterValues[col] = { isNum: true, val: [-Number.MAX_VALUE, Number.MAX_VALUE] };
          } else {
            tempFilterValues[col] = { isNum: false, val: [] };
          }
        }
      });
      for (const i in Object.keys(tempFilterValues)) {
        // remove constraints
        const col = Object.keys(tempFilterValues)[i];
        if (!constraintCols.includes(col)) {
          delete tempFilterValues[col];
        }
      }
      setFilterValues(tempFilterValues);
    }
    // eslint-disable-next-line
  }, [constraintCols, dataset.columns]);

  return (
    <div>
      <Box paddingTop={2} paddingRight={2}>
        {Object.keys(filterValues).map((key) => {
          const value = filterValues[key];

          if (value.isNum) {
            return (
              <RangeFilter
                key={key}
                dataset={dataset}
                col={key}
                value={value.val}
                setValue={(newValue) => {
                  const tempFilterValues = { ...filterValues };
                  tempFilterValues[key].val = newValue;
                  setFilterValues(tempFilterValues);
                }}
                remove={removeFilter}
              />
            );
          }
          return (
            <CategoryFilter
              key={key}
              dataset={dataset}
              col={key}
              value={value.val}
              setValue={(newValue) => {
                const tempFilterValues = { ...filterValues };
                tempFilterValues[key].val = newValue;
                setFilterValues(tempFilterValues);
              }}
              remove={removeFilter}
            />
          );
        })}
      </Box>

      <Box paddingLeft={2} paddingTop={1}>
        <Tooltip
          placement="right"
          title={
            <Typography variant="subtitle2">
              Save the current filter configuration and update the subset of {state.globalLabels.itemLabelPlural}, shown in the web application, accordingly.
            </Typography>
          }
        >
          <Button
            fullWidth
            variant="outlined"
            onClick={() => {
              updateBackendConstraints(filterValues, dataset, triggerDatasetUpdate, state);
            }}
          >
            <FilterAltIcon />
            &nbsp;Apply filter
          </Button>
        </Tooltip>
      </Box>
      <Box paddingLeft={2} paddingTop={1}>
        <Tooltip
          placement="right"
          title={
            <Typography variant="subtitle2">
              Reset filter to the initial configuration. This removes all filters and sets the &quot;experiment cycle&quot; filter to be positive.
            </Typography>
          }
        >
          <Button
            fullWidth
            variant="outlined"
            aria-label="Reset filter to initial state"
            onClick={() => {
              ReactionCIMEBackendFromEnv.resetPOIConstraints(dataset.info.path).then((res_constraints) => {
                if (triggerDatasetUpdate != null) {
                  triggerDatasetUpdate({
                    display: dataset.info.path,
                    path: dataset.info.path,
                    type: dataset.info.type,
                    uploaded: true,
                  });
                }
              });
            }}
          >
            <RotateLeftIcon />
            &nbsp;Reset filter
          </Button>
        </Tooltip>
      </Box>
    </div>
  );
}
