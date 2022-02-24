import { createSlice, PayloadAction } from '@reduxjs/toolkit'


export interface MouseInteractionHooksState {
    mousemove: {x: number, y: number, event_used: boolean}
}

const initialState = { mousemove: null } as MouseInteractionHooksState

const mouseInteractionHooksSlice = createSlice({
    name: 'mouseInteractionHooks',
    initialState,
    reducers:{
        setMouseMove(state, action: PayloadAction<any>) {
            state.mousemove = action.payload
        },
    }
})

export const { setMouseMove } = mouseInteractionHooksSlice.actions
export const mouseInteractionHooks = mouseInteractionHooksSlice.reducer