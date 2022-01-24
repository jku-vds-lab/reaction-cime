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
    const svgRef = React.useRef();
    React.useEffect(() => {
        
        if(legend != null && gRef.current && svgRef.current){
            const rel_width = 250;
            legend.size(rel_width)

            let padding_left, padding_top, padding_right, padding_bottom, rel_height;
            if(legend.height == null){
                rel_height = rel_width;
                padding_left = rel_width*0.1;
                padding_top = padding_left*2;
                padding_right = padding_left*3;
                padding_bottom = padding_left*3;
            }else{
                rel_height = legend.height();
                padding_left = rel_width*0.05;
                padding_top = 0;
                padding_right = padding_left*2;
                padding_bottom = rel_height*1.2;
            }
            
            const gElement = d3.select(gRef.current)
            gElement.html(""); // clear g element
            gElement.call(legend); // draw color legend
            d3.select(svgRef.current).attr("viewBox", `-${padding_left} -${padding_top} ${rel_width+padding_right} ${rel_height+padding_bottom}`)
        }
    }, [legend, gRef, svgRef])
    
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
        <svg ref={svgRef} onClick={()=> {toggleUseVSUP()}} style={{cursor:"pointer"}}><g ref={gRef}></g></svg>
    </>
      
})
  
  