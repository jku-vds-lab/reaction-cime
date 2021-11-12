import { AggregateDataset } from "../Overrides/AggregationTabPanel/AggregateDataset";

const SET = "ducks/aggregatedataset/SET"

interface SetAggregateDatasetAction {
    type: typeof SET
    dataset: AggregateDataset
}


type DatasetActionTypes = SetAggregateDatasetAction

export function setAggregateDatasetAction(dataset: AggregateDataset): DatasetActionTypes {
    return {
        type: SET,
        dataset: dataset
    }
}


const initialState: AggregateDataset = null

export default function aggregateDataset(state = initialState, action): AggregateDataset {
    switch (action.type) {
        case SET:
            return action.dataset
        default:
            return state
    }
}