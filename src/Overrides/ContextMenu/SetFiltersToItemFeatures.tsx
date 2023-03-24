import * as React from 'react';
import { MenuItem } from '@mui/material';
import { EXCLUDED_COLUMNS_ALL, TypedObject, highlightTab } from 'projection-space-explorer';
import { connect, ConnectedProps } from 'react-redux';
import { AppState } from '../../State/Store';
import { updateBackendConstraints } from '../FilterTabPanel/FilterSettings';

const mapStateToProps = (state: AppState) => ({
  globalLabels: state.globalLabels,
  dataset: state.dataset,
  triggerDatasetUpdate: state.handleDataset.triggerUpdate,
  state,
});

const mapDispatchToProps = (dispatch) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  handleClose: () => void;
  pos_x: number;
  pos_y: number;
  menuTarget: TypedObject;
};

export const SetFiltersToItemFeatures = connector(({ handleClose, globalLabels, menuTarget, dataset, triggerDatasetUpdate, state }: Props) => {
  if (menuTarget == null) {
    return null;
  }

  return (
    <MenuItem
      data-cy="set-filters-context-menu"
      onClick={() => {
        const cols = Object.keys(dataset.columns).filter((key) => {
          return dataset.columns[key].project && !EXCLUDED_COLUMNS_ALL.includes(key);
        });
        const constraints = {};
        cols.forEach((col) => {
          const val = [];
          if (dataset.columns[col].isNumeric) {
            const eps = menuTarget[col] * 0.01;
            val.push(menuTarget[col] - eps);
            val.push(menuTarget[col] + eps);
          } else {
            val.push(menuTarget[col]);
          }

          constraints[col] = { val, isNum: dataset.columns[col].isNumeric };
        });

        updateBackendConstraints(constraints, dataset, triggerDatasetUpdate, state);

        highlightTab(7);
        handleClose();
      }}
    >
      {`Set filters to feature values of this ${globalLabels.itemLabel}`}
    </MenuItem>
  );
});
