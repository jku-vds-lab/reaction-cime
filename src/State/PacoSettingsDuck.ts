import { createSlice, PayloadAction } from '@reduxjs/toolkit';
import React from 'react';

export interface PacoSettingsState {
  pacoAttributes: { feature: string; show: boolean; group: string }[];
  pacoConstraints: [];
  pacoRef;
}

const initialState = { pacoAttributes: null, pacoConstraints: [], pacoRef: null } as PacoSettingsState;

const pacoSettingsSlice = createSlice({
  name: 'pacoSettings',
  initialState,
  reducers: {
    setPacoAttributes(state, action: PayloadAction<{ feature: string; show: boolean; group: string }[]>) {
      state.pacoAttributes = action.payload;
    },
    setPacoConstraints(state, action: PayloadAction<[]>) {
      state.pacoConstraints = action.payload;
    },
    setPacoRef(state, action: PayloadAction<React.RefObject<any>>) {
      state.pacoRef = action.payload;
    },
  },
});

export const { setPacoAttributes, setPacoConstraints, setPacoRef } = pacoSettingsSlice.actions;
export const PacoActions = { ...pacoSettingsSlice.actions };
export const pacoSettings = pacoSettingsSlice.reducer;
