import { createSlice, PayloadAction } from '@reduxjs/toolkit'

interface AggregateColorState {
    value_col: string
    uncertainty_col: string
    cache_cols: string[]
}

const initialState = { value_col: "None", uncertainty_col: "None", cache_cols:null } as AggregateColorState

const aggregateColorSlice = createSlice({
    name: 'aggregateColor',
    initialState,
    reducers:{
        setAggregateColor: {
            reducer: (state, action: PayloadAction<AggregateColorState>) => {
                state.value_col = action.payload.value_col
                state.uncertainty_col = action.payload.uncertainty_col
                state.cache_cols = action.payload.cache_cols
            },
            prepare: (values: AggregateColorState) => {
                return {payload: values ?? initialState }
            }
        }
    }
})

export const { setAggregateColor } = aggregateColorSlice.actions
export default aggregateColorSlice.reducer

// const SET = "ducks/aggregateColor/SET"

// const aggregateColor = (state = null, action) => {
//     switch (action.type) {
//         case SET:
//             return action.aggregateColor
//         default:
//             return state
//     }
// }


// export const setAggregateColor = aggregateColor => ({
//     type: SET,
//     aggregateColor: aggregateColor
// })

// export default aggregateColor