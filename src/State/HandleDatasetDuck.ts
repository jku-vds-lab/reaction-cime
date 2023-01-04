import { createSlice, PayloadAction } from '@reduxjs/toolkit';
import { DatasetType } from 'projection-space-explorer';

export interface HandleDatasetState {
  triggerUpdate: (entry: { display: string; path: string; type: DatasetType; uploaded: boolean }, state?) => void;
}

const initialState = { triggerUpdate: null } as HandleDatasetState;

const handleDatasetSlice = createSlice({
  name: 'handleDataset',
  initialState,
  reducers: {
    setTriggerUpdate(state, action: PayloadAction<any>) {
      state.triggerUpdate = action.payload;
    },
  },
});

export const { setTriggerUpdate } = handleDatasetSlice.actions;
export const handleDataset = handleDatasetSlice.reducer;
