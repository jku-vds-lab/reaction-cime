import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export interface SelectionState {
  currentAggregateSelection: any;
}

const initialState = { currentAggregateSelection: null } as SelectionState;

const selectionSlice = createSlice({
  name: 'selection',
  initialState,
  reducers: {
    setCurrentAggregateSelection(state, action: PayloadAction<any>) {
      state.currentAggregateSelection = action.payload;
    },
  },
});

export const { setCurrentAggregateSelection } = selectionSlice.actions;
export const selection = selectionSlice.reducer;
