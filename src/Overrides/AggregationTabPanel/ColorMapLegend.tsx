import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import * as d3 from 'd3v5';
import { InputLabel, MenuItem, Select } from "@mui/material";
import { D3_CONTINUOUS_COLOR_SCALE_LIST, setAggregateColorScale, toggleUseVSUP } from "../../State/AggregateSettingsDuck";

const mapStateToProps = (state: AppState) => ({
    legend: state.aggregateSettings?.legend,
    colorScale: state.aggregateSettings?.colorscale
});

const mapDispatchToProps = (dispatch) => ({
    setAggregateColorScale: value => dispatch(setAggregateColorScale(value)),
    toggleUseVSUP: () => dispatch(toggleUseVSUP())
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    selectAttribute: {key:string, name:string, col_info:any}
};
  

  
export const ColorMapLegend = connector(({legend, colorScale, setAggregateColorScale, selectAttribute, toggleUseVSUP}: Props) => {
    if(selectAttribute == null || selectAttribute.key === "None" || selectAttribute.key == null){
        return null;
    }
    // console.log(Object.keys(d3)) // keys that start with "interpolate"

    const discrete_steps = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

    const handleChange = (event) => {
        setAggregateColorScale(event.target.value);
    };
  
    // add colormap legend
    const gRef = React.useRef();
    React.useEffect(() => {
        if(legend != null && gRef.current){
            const gElement = d3.select(gRef.current)
            gElement.html(""); // clear g element
            gElement.call(legend); // draw color legend
        }
    }, [legend, gRef])
    
    return <>
        <InputLabel id="colorscale-select-label">Choose Colormap</InputLabel>
        <Select
          labelId="colorscale-select-label"
          id="colorscale-select"
          value={colorScale}
          onChange={handleChange}
        >
            {D3_CONTINUOUS_COLOR_SCALE_LIST.map((colorscale) => 
                <MenuItem key={colorscale} value={colorscale} title={colorscale}>
                    <div style={{ width: '100%', minWidth: '13rem', height: '1rem', backgroundImage: `linear-gradient(to right, ${discrete_steps.map(step => d3[colorscale](step)).join(',')})` }}></div>
                </MenuItem>)
            }
        </Select>
        <svg onClick={()=> {toggleUseVSUP()}} style={{width:"100%", height:"300px", cursor:"pointer"}}><g ref={gRef}></g></svg>
    </>
      
})
  
  