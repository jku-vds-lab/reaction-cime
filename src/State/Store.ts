import { createAction } from "@reduxjs/toolkit";
import { createRootReducer, createViewDuckReducer, RootActionTypes, RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import aggregateSettings, { AggregateInitStates } from "./AggregateSettingsDuck";
import { handleDataset } from "./HandleDatasetDuck";
// import Dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";
// import cimeBackgroundSelection from "projection-space-explorer/dist/components/Ducks/CimeBackgroundSelectionDuck";
import lineUpInput from "./LineUpInputDuck";
import { mouseInteractionHooks } from "./MouseInteractionHooksDuck";
import { pacoSettings } from "./PacoSettingsDuck";
import { selection } from "./SelectionDuck";

const resetViews = createAction<void>('view/resetView');

export const CIME4RViewActions = {
  resetViews
}

export const CIMEReducers = {
  lineUpInput: lineUpInput,
  // cimeBackgroundSelection: cimeBackgroundSelection,
  // dataset: Dataset,
  // aggregateSettings: aggregateSettings,
  multiples: createViewDuckReducer({ aggregateSettings }, (builder) => {
    builder.addCase(resetViews, (state, action) => {
      for(let i in state.multiples.ids){
        const id = state.multiples.ids[i]
        state.multiples.entities[id].attributes.aggregateSettings = AggregateInitStates.all
      }
    })
  }).reducer,
  mouseInteractionHooks: mouseInteractionHooks,
  selection: selection,
  handleDataset: handleDataset,
  pacoSettings: pacoSettings
};

const cimeCombined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof cimeCombined>;

export type AppState = RootState & CimeState;

export function createCIMERootReducer() {
  const pseRootReducer = createRootReducer(CIMEReducers)
  
  return (state: Parameters<typeof pseRootReducer>[0], action: Parameters<typeof pseRootReducer>[1]) => {
    const newState = pseRootReducer(state, action)
    // console.log(action.type)
    if (action.type === RootActionTypes.DATASET) {
      // initialize pacoAttributes when dataset changes
      // TODO: might not be necessary
      if(newState.dataset != null){
        const newPacoAttributes = Object.keys(newState.dataset.columns).map((col) => {
          return {feature: col, show: newState.dataset.columns[col].metaInformation.paco === true};
        })
        Object.assign(newState, {pacoSettings: {...state.pacoSettings, pacoAttributes: newPacoAttributes}});
      }
    }

    return newState;
  };
}