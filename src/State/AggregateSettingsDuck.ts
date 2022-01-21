import { createSlice, PayloadAction } from '@reduxjs/toolkit'

export interface AggregateSettingsState {
    legend: any
}

const initialState = { legend: null } as AggregateSettingsState

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
        }
    }
})

export const { setAggregateColorMapLegend } = aggregateSettingsSlice.actions
export default aggregateSettingsSlice.reducer
