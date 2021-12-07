import { RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import aggregateDataset from "./AggregateDatasetDuck";
import aggregateColor from "./AggregateColorDuck";
import lineUpInput from "./LineUpInputDuck";

export const CIMEReducers = {
  lineUpInput: lineUpInput,
  aggregateDataset: aggregateDataset,
  aggregateColor: aggregateColor
};

const combined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;

export type AppState = RootState & CimeState;
