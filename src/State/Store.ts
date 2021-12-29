import { RootState } from "projection-space-explorer";
import { combineReducers } from "redux";
import aggregateDataset from "./AggregateDatasetDuck";
import aggregateColor from "./AggregateColorDuck";
import Dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";
import lineUpInput from "./LineUpInputDuck";

export const CIMEReducers = {
  lineUpInput: lineUpInput,
  aggregateDataset: aggregateDataset,
  aggregateColor: aggregateColor,
  dataset: Dataset
};

const combined = combineReducers(CIMEReducers);

/**
 * Cime typings...
 */
export type CimeState = ReturnType<typeof combined>;

export type AppState = RootState & CimeState;
