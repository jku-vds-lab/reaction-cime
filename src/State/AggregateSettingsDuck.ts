import { createSlice, PayloadAction } from '@reduxjs/toolkit'

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
}

const initialState = { legend: null, colorscale: D3_CONTINUOUS_COLOR_SCALE_LIST[0], useVSUP: true, sampleSize: 200 } as AggregateSettingsState

const aggregateSettingsSlice = createSlice({
    name: 'aggregateSettings',
    initialState,
    reducers:{
        setAggregateColorMapLegend(state, action: PayloadAction<any>) {
            state.legend = action.payload
        },
        setAggregateColorScale(state, action: PayloadAction<string>) {
            state.colorscale = action.payload
        },
        toggleUseVSUP(state){
            state.useVSUP = !state.useVSUP
        },
        setSampleSize(state, action: PayloadAction<number>){
            state.sampleSize = action.payload
        }
    }
})

export const { setAggregateColorMapLegend, setAggregateColorScale, toggleUseVSUP, setSampleSize } = aggregateSettingsSlice.actions
export default aggregateSettingsSlice.reducer
