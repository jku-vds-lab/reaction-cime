import { createSlice, PayloadAction } from '@reduxjs/toolkit'

// console.log(Object.keys(d3)) // could also use keys of d3.js that start with "interpolate"
export const D3_CONTINUOUS_COLOR_SCALE_LIST = [
    "interpolateViridis",
    "interpolateCividis",
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


export interface AggregateSettingsState {
    legend: any
    colorscale: string
    useVSUP: boolean
    sampleSize: number
    valueFilter: string[] // contains colorvalue of which values to show; if empty, all are shown;
}

const initialState = { legend: null, colorscale: D3_CONTINUOUS_COLOR_SCALE_LIST[0], useVSUP: true, sampleSize: 200, valueFilter: [] } as AggregateSettingsState

const aggregateSettingsSlice = createSlice({
    name: 'aggregateSettings',
    initialState,
    reducers:{
        setAggregateColorMapLegend(state, action: PayloadAction<any>) {
            state.valueFilter = []
            state.legend = action.payload
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
        }
    }
})

export const { setAggregateColorMapLegend, setAggregateColorScale, toggleUseVSUP, setSampleSize, addValueFilter, removeValueFilter, clearValueFilter } = aggregateSettingsSlice.actions
export default aggregateSettingsSlice.reducer
