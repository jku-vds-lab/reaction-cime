import { createAction } from '@reduxjs/toolkit';
import { createRootReducer, createViewDuckReducer, RootActionTypes, RootState } from 'projection-space-explorer';
import { combineReducers } from 'redux';
import aggregateSettings, { AggregateInitStates } from './AggregateSettingsDuck';
import { handleDataset } from './HandleDatasetDuck';
// import Dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";
// import cimeBackgroundSelection from "projection-space-explorer/dist/components/Ducks/CimeBackgroundSelectionDuck";
import lineUpInput from './LineUpInputDuck';
import { mouseInteractionHooks } from './MouseInteractionHooksDuck';
import { pacoSettings } from './PacoSettingsDuck';
import { projects } from './ProjectsDuck';
import { selection } from './SelectionDuck';

const resetViews = createAction<void>('view/resetView');

export const CIME4RViewActions = {
  resetViews,
};

export const CIMEReducers = {
  lineUpInput,
  // cimeBackgroundSelection: cimeBackgroundSelection,
  // dataset: Dataset,
  // aggregateSettings: aggregateSettings,
  multiples: createViewDuckReducer({ aggregateSettings }, (builder) => {
    builder.addCase(resetViews, (state, action) => {
      state.multiples.ids.forEach((id) => {
        state.multiples.entities[id].attributes.aggregateSettings = AggregateInitStates.all;
      });
    });
  }).reducer,
  mouseInteractionHooks,
  selection,
  handleDataset,
  pacoSettings,
  projects: projects.reducer,
};

const cimeCombined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof cimeCombined>;

export type AppState = RootState & CimeState;

export function createCIMERootReducer() {
  const pseRootReducer = createRootReducer(CIMEReducers);

  return (state: Parameters<typeof pseRootReducer>[0], action: Parameters<typeof pseRootReducer>[1]) => {
    const newState = pseRootReducer(state, action);
    // console.log(action.type)
    if (action.type === RootActionTypes.DATASET) {
      // initialize pacoAttributes when dataset changes
      // TODO: might not be necessary
      if (newState.dataset != null) {
        const newPacoAttributes = Object.keys(newState.dataset.columns).map((col) => {
          return { feature: col, show: newState.dataset.columns[col].metaInformation.paco as boolean };
        });
        Object.assign(newState, { pacoSettings: { ...state.pacoSettings, pacoAttributes: newPacoAttributes } });
      }
    }

    return newState;
  };
}
