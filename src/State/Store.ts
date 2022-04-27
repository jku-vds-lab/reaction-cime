import { createAction } from "@reduxjs/toolkit";
import { createViewDuckReducer, RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import aggregateSettings, { AggregateInitStates } from "./AggregateSettingsDuck";
import { handleDataset } from "./HandleDatasetDuck";
// import Dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";
// import cimeBackgroundSelection from "projection-space-explorer/dist/components/Ducks/CimeBackgroundSelectionDuck";
import lineUpInput from "./LineUpInputDuck";
import { mouseInteractionHooks } from "./MouseInteractionHooksDuck";
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
};

const combined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;

export type AppState = RootState & CimeState;
