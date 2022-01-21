import { Slider, Typography } from "@mui/material";
import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import { AppState } from "../../State/Store";

const mapStateToPropsSlider = (state: AppState) => ({
    workspace: state.projections.workspace,
    poiDataset: state.dataset
  });
  
  const mapDispatchToPropsSlider = (dispatch) => ({
    setAggregateColor: values => dispatch(setAggregateColor(values)),
  });
  
  const sliderconnector = connect(mapStateToPropsSlider, mapDispatchToPropsSlider);
  
  type SliderPropsFromRedux = ConnectedProps<typeof sliderconnector>;
  
  type SliderProps = SliderPropsFromRedux & {
    selectAttribute: {key:string, name:string, col_info:any}
  };
  
  
  
export const StepSlider = sliderconnector(({selectAttribute, setAggregateColor, workspace, poiDataset}: SliderProps) => {
    if(selectAttribute == null || selectAttribute.key === "None" || selectAttribute.key == null){
        return null;
    }

    React.useEffect(() => {
        // reset aggregate color to hide the aggregated dataset in the background
        setAggregateColor(null)
    // eslint-disable-next-line
    }, [workspace, poiDataset]) // this is triggered during the embedding


    const stepChanged = (newVal) => {
        setCurStep(newVal)
    }

    const [curStep, setCurStep] = React.useState(0);
    const [marks, setMarks] = React.useState([]);


    React.useEffect(() => {
        let m;
        if(selectAttribute.col_info){
        const variables = Object.keys(selectAttribute.col_info);
        let step_arr = Object.keys(selectAttribute.col_info[variables[0]].temporal_columns)
        step_arr.sort()

        m = step_arr.map((step, index) => {return {value: parseInt(step), label: step}})
        }
        setMarks(m)
        // setCurStep(m[m.length-1].value)
    }, [selectAttribute])

    React.useEffect(() => {
        const timestep = curStep;
        if(selectAttribute.col_info){
        const variables_array = Object.keys(selectAttribute.col_info);
        // TODO: users should be able to select which column is value and which is uncertainty
        let value_col = selectAttribute.col_info[variables_array[0]]["temporal_columns"][timestep]
        let cache_cols = Object.values(selectAttribute.col_info[variables_array[0]]["temporal_columns"])

        let uncertainty_col = null;
        if(variables_array.length >= 2){
            uncertainty_col = selectAttribute.col_info[variables_array[1]]["temporal_columns"][timestep]
            cache_cols = cache_cols.concat(Object.values(selectAttribute.col_info[variables_array[1]]["temporal_columns"]))
        }
        
        setAggregateColor({"value_col": value_col, "uncertainty_col": uncertainty_col, "cache_cols": cache_cols});
        }else{
        setAggregateColor({"value_col": selectAttribute.key, "uncertainty_col": null, "cache_cols": null});
        }
        // eslint-disable-next-line
    }, [curStep, selectAttribute])

    return (marks && marks.length > 0) && <>
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
    </>
})
  
  