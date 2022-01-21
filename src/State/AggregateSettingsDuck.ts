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
}

const initialState = { legend: null, colorscale: D3_CONTINUOUS_COLOR_SCALE_LIST[0], useVSUP: true } as AggregateSettingsState

const aggregateSettingsSlice = createSlice({
    name: 'aggregateSettings',
    initialState,
    reducers:{
        setAggregateColorMapLegend: {
            reducer: (state, action: PayloadAction<any>) => {
                state.legend = action.payload
            },
            prepare: (value: any) => {
                return {payload: value }
            }
        },
        setAggregateColorScale: {
            reducer: (state, action: PayloadAction<string>) => {
                state.colorscale = action.payload
            },
            prepare: (value: any) => {
                return {payload: value }
            }
        },
        toggleUseVSUP(state){
            state.useVSUP = !state.useVSUP
        }
    }
})

export const { setAggregateColorMapLegend, setAggregateColorScale, toggleUseVSUP } = aggregateSettingsSlice.actions
export default aggregateSettingsSlice.reducer
