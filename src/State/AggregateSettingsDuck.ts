import { createSlice, PayloadAction } from '@reduxjs/toolkit';

export enum AggregationMethod {
  MIN = 'min',
  MAX = 'max',
  MEDIAN = 'median',
  MEAN = 'mean',
}

// console.log(Object.keys(d3)) // could also use keys of d3.js that start with "interpolate"
export const D3_CONTINUOUS_COLOR_SCALE_LIST = [
  'interpolateCividis',
  'interpolateViridis',
  'interpolateYlGnBu',
  'interpolateYlOrBr',
  'interpolateMagma',
  'interpolateInferno',
  'interpolateWarm',
  'interpolateCool',
  'interpolateBuGn',
  'interpolateGnBu',
  'interpolateOrRd',
  'interpolateYlOrRd',
  'interpolateTurbo',
];

type AggregateColorType = {
  value_col: string;
  uncertainty_col: string;
  cache_cols: string[];
};
const initialAggregateColor = { value_col: 'None', uncertainty_col: 'None', cache_cols: null } as AggregateColorType;

type ColorMapSettingsType = {
  aggregateColor: AggregateColorType;
  scale_obj: any;
  colorscale: string;
  useVSUP: boolean;
  valueFilter: string[]; // contains colorvalue of which values to show; if empty, all are shown;
};
const initialColorMapSettings = {
  aggregateColor: initialAggregateColor,
  scale_obj: null,
  colorscale: D3_CONTINUOUS_COLOR_SCALE_LIST[0],
  useVSUP: true,
  valueFilter: [],
} as ColorMapSettingsType;

type StepperSettingsType = {
  curStep: number;
};
const initialStepperSettings = { curStep: 0 } as StepperSettingsType;

type AdvancedSettingsType = {
  deriveRange: boolean;
  valueRange: { min: number; max: number };
  uncertaintyRange: { min: number; max: number };
  variableIndex: { valueVariableIndex: number; uncertaintyVariableIndex: number };
  aggregationMethod: { valueAggregationMethod: string; uncertaintyAggregationMethod: string };
};
const initialAdvancedSettings = {
  deriveRange: true,
  valueRange: null,
  uncertaintyRange: null,
  variableIndex: { valueVariableIndex: 0, uncertaintyVariableIndex: 0 },
  aggregationMethod: { valueAggregationMethod: AggregationMethod.MAX, uncertaintyAggregationMethod: AggregationMethod.MIN },
  sampleSize: 20,
} as AdvancedSettingsType;

export interface AggregateSettingsState {
  colormapSettings: ColorMapSettingsType;
  stepperSettings: StepperSettingsType;
  advancedSettings: AdvancedSettingsType;

  selectAttribute: string;
  selectAttributeInfo: {};
}

const initialState = {
  colormapSettings: initialColorMapSettings,
  stepperSettings: initialStepperSettings,
  advancedSettings: initialAdvancedSettings,
  selectAttribute: null,
  selectAttributeInfo: null,
} as AggregateSettingsState;

export const AggregateInitStates = {
  all: initialState,
  colormapSettings: initialColorMapSettings,
  stepperSettings: initialStepperSettings,
  advancedSettings: initialAdvancedSettings,
};

const aggregateSettingsSlice = createSlice({
  name: 'aggregateSettings',
  initialState,
  reducers: {
    setSelectAttribute(state, action: PayloadAction<{ attribute_name: string; attribute_info: {} }>) {
      state.selectAttribute = action.payload.attribute_name;
      state.selectAttributeInfo = action.payload.attribute_info;
    },
    setCurStep(state, action: PayloadAction<number>) {
      state.stepperSettings.curStep = action.payload;
    },
    setAggregateColor: {
      reducer: (state, action: PayloadAction<AggregateColorType>) => {
        if (
          state.colormapSettings.aggregateColor.value_col === action.payload.value_col &&
          state.colormapSettings.aggregateColor.uncertainty_col === action.payload.uncertainty_col
        ) {
          return;
        }
        state.colormapSettings.scale_obj = null;
        state.colormapSettings.valueFilter = [];

        state.colormapSettings.aggregateColor.value_col = action.payload.value_col;
        state.colormapSettings.aggregateColor.uncertainty_col = action.payload.uncertainty_col;
        state.colormapSettings.aggregateColor.cache_cols = action.payload.cache_cols;
      },
      prepare: (values: AggregateColorType) => {
        return { payload: values ?? initialAggregateColor };
      },
    },
    setAggregateColorMapScale(state, action: PayloadAction<any>) {
      state.colormapSettings.valueFilter = [];
      state.colormapSettings.scale_obj = action.payload;
    },
    setValueRange(state, action: PayloadAction<{ min: number; max: number }>) {
      state.advancedSettings.valueRange = action.payload;
    },
    setUncertaintyRange(state, action: PayloadAction<{ min: number; max: number }>) {
      state.advancedSettings.uncertaintyRange = action.payload;
    },
    toggleDeriveRange(state) {
      state.advancedSettings.deriveRange = !state.advancedSettings.deriveRange;
    },
    setDeriveRange(state, action: PayloadAction<boolean>) {
      state.advancedSettings.deriveRange = action.payload;
    },
    setAggregateColorScale(state, action: PayloadAction<string>) {
      state.colormapSettings.valueFilter = [];
      state.colormapSettings.colorscale = action.payload;
    },
    toggleUseVSUP(state) {
      state.colormapSettings.valueFilter = [];
      state.colormapSettings.useVSUP = !state.colormapSettings.useVSUP;
    },
    addValueFilter(state, action: PayloadAction<string>) {
      state.colormapSettings.valueFilter.push(action.payload);
    },
    removeValueFilter(state, action: PayloadAction<string>) {
      const index = state.colormapSettings.valueFilter.indexOf(action.payload);
      if (index !== -1) {
        state.colormapSettings.valueFilter.splice(index, 1);
      }
    },
    clearValueFilter(state) {
      state.colormapSettings.valueFilter = [];
    },
    setVariableIndex(state, action: PayloadAction<{ valueVariableIndex: number; uncertaintyVariableIndex: number }>) {
      state.advancedSettings.variableIndex = action.payload;
    },
    setAggregationMethod(state, action: PayloadAction<{ valueAggregationMethod: string; uncertaintyAggregationMethod: string }>) {
      state.advancedSettings.aggregationMethod = action.payload;
    },
  },
});

export const {
  setCurStep,
  setSelectAttribute,
  setAggregateColor,
  setAggregateColorMapScale,
  setValueRange,
  setUncertaintyRange,
  toggleDeriveRange,
  setDeriveRange,
  setAggregateColorScale,
  toggleUseVSUP,
  addValueFilter,
  removeValueFilter,
  clearValueFilter,
  setVariableIndex,
  setAggregationMethod,
} = aggregateSettingsSlice.actions;
export default aggregateSettingsSlice.reducer;
