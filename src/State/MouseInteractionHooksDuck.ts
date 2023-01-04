import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export interface MouseInteractionHooksState {
  mousemove: { x: number; y: number; event_used: boolean };
  mousedown: { x: number; y: number; timestamp: number };
  mouseclick: { x: number; y: number; event_used: boolean; button: number };
}

const initialState = { mousemove: null, mousedown: null, mouseclick: null } as MouseInteractionHooksState;

const mouseInteractionHooksSlice = createSlice({
  name: 'mouseInteractionHooks',
  initialState,
  reducers: {
    setMouseMove(state, action: PayloadAction<{ x: number; y: number; event_used: boolean }>) {
      state.mousemove = action.payload;
    },
    // setMouseDown(state, action: PayloadAction<{x: number, y: number}>){
    //     state.mousedown = {x: action.payload.x, y: action.payload.y, timestamp: Date.now()}
    //     state.mouseclick = null
    // },
    // setMouseUp(state, action: PayloadAction<{x: number, y: number}>){
    //     if(state?.mousedown?.timestamp != null){
    //         const delta_time = Date.now() - state.mousedown.timestamp
    //         // if time between mousedown and mouseup is smaller than 1 sec and the positions of mousedown and mouseup are the same, we count it as click
    //         if(delta_time < 1000 && state.mousedown.x === action.payload.x && state.mousedown.y === action.payload.y){
    //             state.mousedown = null
    //             state.mouseclick = {...action.payload, event_used=false}
    //         }
    //     }

    // },
    setMouseClick(state, action: PayloadAction<{ x: number; y: number; event_used: boolean; button: number }>) {
      state.mouseclick = action.payload;
    },
  },
});

export const { setMouseMove, setMouseClick } = mouseInteractionHooksSlice.actions;
export const mouseInteractionHooks = mouseInteractionHooksSlice.reducer;
