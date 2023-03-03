import { createSlice, PayloadAction } from '@reduxjs/toolkit';
import { GenericFingerprintAttribute } from 'projection-space-explorer';
import React from 'react';

export interface PacoSettingsState {
  pacoAttributes: GenericFingerprintAttribute[];
  pacoConstraints: [];
  pacoRef;
}

const initialState = { pacoAttributes: null, pacoConstraints: [], pacoRef: null } as PacoSettingsState;

const pacoSettingsSlice = createSlice({
  name: 'pacoSettings',
  initialState,
  reducers: {
    setPacoAttributes(state, action: PayloadAction<GenericFingerprintAttribute[]>) {
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
