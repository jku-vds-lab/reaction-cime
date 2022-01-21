import { RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import aggregateColor from "./AggregateColorDuck";
import aggregateSettings from "./AggregateSettingsDuck";
import Dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";
import cimeBackgroundSelection from "projection-space-explorer/dist/components/Ducks/CimeBackgroundSelectionDuck";
import lineUpInput from "./LineUpInputDuck";

export const CIMEReducers = {
  lineUpInput: lineUpInput,
  aggregateColor: aggregateColor,
  cimeBackgroundSelection: cimeBackgroundSelection,
  dataset: Dataset,
  aggregateSettings: aggregateSettings
};

const combined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;

export type AppState = RootState & CimeState;
