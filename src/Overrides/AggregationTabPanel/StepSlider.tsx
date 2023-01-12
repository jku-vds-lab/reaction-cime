import { Slider, Typography } from '@mui/material';
import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { AggregateActions } from '../../State/AggregateSettingsDuck';
import { AppState } from '../../State/Store';

const mapStateToPropsSlider = (state: AppState) => ({
  variableIndex: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.advancedSettings.variableIndex,
  selectAttribute: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttribute,
  selectAttributeInfo: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttributeInfo,
  curStep: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.stepperSettings.curStep,
});

const mapDispatchToPropsSlider = (dispatch) => ({
  setAggregateColor: (values) => dispatch(AggregateActions.setAggregateColor(values)),
  setCurStep: (values) => dispatch(AggregateActions.setCurStep(values)),
});

const sliderconnector = connect(mapStateToPropsSlider, mapDispatchToPropsSlider);

type SliderPropsFromRedux = ConnectedProps<typeof sliderconnector>;

type SliderProps = SliderPropsFromRedux;

export const StepSlider = sliderconnector(({ selectAttribute, setAggregateColor, selectAttributeInfo, variableIndex, curStep, setCurStep }: SliderProps) => {
  const stepChanged = (newVal) => {
    setCurStep(newVal);
  };

  const [marks, setMarks] = React.useState([]);

  React.useEffect(() => {
    let m;
    if (selectAttributeInfo) {
      const variables = Object.keys(selectAttributeInfo);
      const stepArr = Object.keys(selectAttributeInfo[variables[0]].temporal_columns);
      stepArr.sort();

      m = stepArr.map((step, index) => {
        return { value: parseInt(step, 10), label: step };
      });
    }
    setMarks(m);
    // setCurStep(m[m.length-1].value)
  }, [selectAttribute, selectAttributeInfo]);

  React.useEffect(() => {
    const timestep = curStep;
    if (selectAttributeInfo && variableIndex) {
      const variablesArray = Object.keys(selectAttributeInfo);
      const valueCol = selectAttributeInfo[variablesArray[variableIndex.valueVariableIndex]].temporal_columns[timestep];
      let cacheCols = Object.values(selectAttributeInfo[variablesArray[variableIndex.valueVariableIndex]].temporal_columns);

      let uncertaintyCol = null;
      if (variableIndex.valueVariableIndex !== variableIndex.uncertaintyVariableIndex && variablesArray[variableIndex.uncertaintyVariableIndex] != null) {
        uncertaintyCol = selectAttributeInfo[variablesArray[variableIndex.uncertaintyVariableIndex]].temporal_columns[timestep];
        cacheCols = cacheCols.concat(Object.values(selectAttributeInfo[variablesArray[variableIndex.uncertaintyVariableIndex]].temporal_columns));
      }
      setAggregateColor({ value_col: valueCol, uncertainty_col: uncertaintyCol, cache_cols: cacheCols });
    } else {
      setAggregateColor({ value_col: selectAttribute, uncertainty_col: null, cache_cols: null });
    }
    // eslint-disable-next-line
  }, [curStep, selectAttribute, variableIndex]);

  return marks && marks.length > 0 ? (
    <>
      <Typography id="range-slider" gutterBottom>
        Choose step
      </Typography>
      <Slider
        min={0}
        max={marks.length - 1}
        value={curStep}
        onChange={(_, newValue) => stepChanged(newValue)}
        step={1}
        marks={marks}
        valueLabelDisplay="auto"
      />
    </>
  ) : null;
});
