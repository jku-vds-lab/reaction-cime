import { Slider, Typography } from "@mui/material";
import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { setAggregateColor } from "../../State/AggregateSettingsDuck";
import { AppState } from "../../State/Store";

const mapStateToPropsSlider = (state: AppState) => ({
    variableIndex: state.aggregateSettings?.variableIndex,
  });
  
  const mapDispatchToPropsSlider = (dispatch) => ({
    setAggregateColor: values => dispatch(setAggregateColor(values)),
  });
  
  const sliderconnector = connect(mapStateToPropsSlider, mapDispatchToPropsSlider);
  
  type SliderPropsFromRedux = ConnectedProps<typeof sliderconnector>;
  
  type SliderProps = SliderPropsFromRedux & {
    selectAttribute: {key:string, name:string}
    selectAttributeInfo: any
  };
  
  
  
export const StepSlider = sliderconnector(({selectAttribute, setAggregateColor, selectAttributeInfo, variableIndex}: SliderProps) => {


    const stepChanged = (newVal) => {
        setCurStep(newVal)
    }

    const [curStep, setCurStep] = React.useState(0);
    const [marks, setMarks] = React.useState([]);


    React.useEffect(() => {
        let m;
        if(selectAttributeInfo){
            const variables = Object.keys(selectAttributeInfo);
            let step_arr = Object.keys(selectAttributeInfo[variables[0]].temporal_columns)
            step_arr.sort()

            m = step_arr.map((step, index) => {return {value: parseInt(step), label: step}})
        }
        setMarks(m)
        // setCurStep(m[m.length-1].value)
    }, [selectAttribute, selectAttributeInfo])

    React.useEffect(() => {
        const timestep = curStep;
        if(selectAttributeInfo && variableIndex){
            const variables_array = Object.keys(selectAttributeInfo);
            let value_col = selectAttributeInfo[variables_array[variableIndex.valueVariableIndex]]["temporal_columns"][timestep]
            let cache_cols = Object.values(selectAttributeInfo[variables_array[variableIndex.valueVariableIndex]]["temporal_columns"])

            let uncertainty_col = null;
            if(variableIndex.valueVariableIndex !== variableIndex.uncertaintyVariableIndex && variables_array[variableIndex.uncertaintyVariableIndex] != null){
                uncertainty_col = selectAttributeInfo[variables_array[variableIndex.uncertaintyVariableIndex]]["temporal_columns"][timestep]
                cache_cols = cache_cols.concat(Object.values(selectAttributeInfo[variables_array[variableIndex.uncertaintyVariableIndex]]["temporal_columns"]))
            }
        
            setAggregateColor({"value_col": value_col, "uncertainty_col": uncertainty_col, "cache_cols": cache_cols});
        }else{
            setAggregateColor({"value_col": selectAttribute.key, "uncertainty_col": null, "cache_cols": null});
        }
        // eslint-disable-next-line
    }, [curStep, selectAttribute, variableIndex])

    return <>{(marks && marks.length > 0) && <>
        <Typography id="range-slider" gutterBottom>
            Choose Step
        </Typography>
        <Slider
            min={0}
            max={marks.length-1}
            value={curStep}
            onChange={(_, newValue) => stepChanged(newValue)}
            step={1}
            marks={marks}
            valueLabelDisplay="auto"
        ></Slider>
    </>}
    </>
})
  
  