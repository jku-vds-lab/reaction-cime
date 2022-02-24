import { createSlice, PayloadAction } from '@reduxjs/toolkit'

export enum AggregationMethod {
    MIN = "min",
    MAX = "max",
    MEDIAN = "median",
    MEAN = "mean"
}


// console.log(Object.keys(d3)) // could also use keys of d3.js that start with "interpolate"
export const D3_CONTINUOUS_COLOR_SCALE_LIST = [
    "interpolateCividis",
    "interpolateViridis",
    "interpolateYlGnBu",
    "interpolateYlOrBr",
    "interpolateMagma",
    "interpolateInferno",
    "interpolateWarm",
    "interpolateCool",
    "interpolateBuGn",
    "interpolateGnBu",
    "interpolateOrRd",
    "interpolateYlOrRd",
    "interpolateTurbo",
]

type AggregateColorType = {
    value_col: string,
    uncertainty_col: string,
    cache_cols: string[]
}

const initialAggregateColor = { value_col: "None", uncertainty_col: "None", cache_cols: null } as AggregateColorType


export interface AggregateSettingsState {
    aggregateColor: AggregateColorType
    scale_obj: any
    deriveRange: boolean
    valueRange: {min: number, max: number}
    uncertaintyRange: {min: number, max: number}
    variableIndex: {valueVariableIndex: number, uncertaintyVariableIndex: number}
    aggregationMethod: {valueAggregationMethod: string, uncertaintyAggregationMethod: string}
    colorscale: string
    useVSUP: boolean
    sampleSize: number
    valueFilter: string[] // contains colorvalue of which values to show; if empty, all are shown;
}

const initialState = { aggregateColor: initialAggregateColor, scale_obj: null, valueRange: null, uncertaintyRange: null, deriveRange: true, colorscale: D3_CONTINUOUS_COLOR_SCALE_LIST[0], useVSUP: true, sampleSize: 200, valueFilter: [], variableIndex: {valueVariableIndex: 0, uncertaintyVariableIndex: 0}, aggregationMethod: {valueAggregationMethod: AggregationMethod.MAX, uncertaintyAggregationMethod: AggregationMethod.MIN} } as AggregateSettingsState

const aggregateSettingsSlice = createSlice({
    name: 'aggregateSettings',
    initialState,
    reducers:{
        setAggregateColor: {
            reducer: (state, action: PayloadAction<AggregateColorType>) => {
                state.scale_obj = null;
                state.valueFilter = [];

                state.aggregateColor.value_col = action.payload.value_col
                state.aggregateColor.uncertainty_col = action.payload.uncertainty_col
                state.aggregateColor.cache_cols = action.payload.cache_cols
            },
            prepare: (values: AggregateColorType) => {
                return {payload: values ?? initialAggregateColor }
            }
        },
        setAggregateColorMapScale(state, action: PayloadAction<any>) {
            state.valueFilter = []
            state.scale_obj = action.payload
        },
        setValueRange(state, action:PayloadAction<{min: number, max: number}>){
            state.valueRange = action.payload
        },
        setUncertaintyRange(state, action:PayloadAction<{min: number, max: number}>){
            state.uncertaintyRange = action.payload
        },
        toggleDeriveRange(state){
            state.deriveRange = !state.deriveRange
        },
        setDeriveRange(state, action: PayloadAction<boolean>){
            state.deriveRange = action.payload
        },
        setAggregateColorScale(state, action: PayloadAction<string>) {
            state.valueFilter = []
            state.colorscale = action.payload
        },
        toggleUseVSUP(state){
            state.valueFilter = []
            state.useVSUP = !state.useVSUP
        },
        setSampleSize(state, action: PayloadAction<number>){
            state.sampleSize = action.payload
        },
        addValueFilter(state, action: PayloadAction<string>){
            state.valueFilter.push(action.payload)
        },
        removeValueFilter(state, action: PayloadAction<string>){
            var index = state.valueFilter.indexOf(action.payload);
            if (index !== -1) {
                state.valueFilter.splice(index, 1);
            }
        },
        clearValueFilter(state){
            state.valueFilter = []
        },
        setVariableIndex(state, action: PayloadAction<{valueVariableIndex: number, uncertaintyVariableIndex: number}>){
            state.variableIndex = action.payload
        },
        setAggregationMethod(state, action: PayloadAction<{valueAggregationMethod: string, uncertaintyAggregationMethod: string}>){
            state.aggregationMethod = action.payload
        },
    }
})

export const { setAggregateColor, setAggregateColorMapScale, setValueRange, setUncertaintyRange, 
    toggleDeriveRange, setDeriveRange, setAggregateColorScale, toggleUseVSUP, setSampleSize, addValueFilter, 
    removeValueFilter, clearValueFilter, setVariableIndex, setAggregationMethod } = aggregateSettingsSlice.actions
export default aggregateSettingsSlice.reducer
