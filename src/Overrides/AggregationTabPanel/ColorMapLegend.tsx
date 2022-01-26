import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import * as d3 from 'd3v5';
import { Button, InputLabel, MenuItem, Select } from "@mui/material";
import { D3_CONTINUOUS_COLOR_SCALE_LIST, setAggregateColorScale, addValueFilter, removeValueFilter, toggleUseVSUP, clearValueFilter } from "../../State/AggregateSettingsDuck";

const mapStateToProps = (state: AppState) => ({
    legend: state.aggregateSettings?.legend,
    colorScale: state.aggregateSettings?.colorscale
});

const mapDispatchToProps = (dispatch) => ({
    setAggregateColorScale: value => dispatch(setAggregateColorScale(value)),
    toggleUseVSUP: () => dispatch(toggleUseVSUP()),
    addValueFilter: (value) => dispatch(addValueFilter(value)),
    removeValueFilter: (value) => dispatch(removeValueFilter(value)),
    clearValueFilter: () => dispatch(clearValueFilter()),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    selectAttribute: {key:string, name:string, col_info:any}
};
  

  
export const ColorMapLegend = connector(({legend, colorScale, setAggregateColorScale, selectAttribute, toggleUseVSUP, addValueFilter, removeValueFilter, clearValueFilter }: Props) => {
    if(selectAttribute == null || selectAttribute.key === "None" || selectAttribute.key == null){
        return null;
    }

    const discrete_steps = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

    const handleChange = (event) => {
        setAggregateColorScale(event.target.value);
    };
  
    // add colormap legend
    const gRef = React.useRef();
    const svgRef = React.useRef();
    const [colorSections, setColorSections] = React.useState([]);

    const clearFilter = (color_sections) => {
        clearValueFilter()
        
        color_sections.attr("stroke", null)
        color_sections.attr("stroke-width", "0px")
    }

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

            const svgElement = d3.select(svgRef.current)
            svgElement.attr("viewBox", `-${padding_left} -${padding_top} ${rel_width+padding_right} ${rel_height+padding_bottom}`)

            // --- add interaction with legend
            const legend_container = svgElement.select(".legend > g:last-child")
            let color_sections = legend_container.selectAll("path")
            if(color_sections.nodes().length <= 0){ // if there are no path elements, we look for rect elements
                color_sections = gElement.selectAll("rect")
            }
            setColorSections(color_sections)
            
            color_sections
            // .on("mouseover", (d, i) => {
            //      //stroke="black" stroke-width="0.02px"
            //     d3.select(color_sections.nodes()[i]).attr("stroke", "black")
            // })
            // .on("mouseout", (d, i) => {
            //     d3.select(color_sections.nodes()[i]).attr("stroke", "none")
            // })
            .on("click", (d, i) => {
                // color_sections.attr("stroke", "none")
                const cur_node = d3.select(color_sections.nodes()[i])
                if(cur_node.attr("stroke") == null){
                    cur_node.attr("stroke", "#1f77b4")
                    cur_node.attr("stroke-width", "3px")
                    addValueFilter(cur_node.attr("fill"))
                }else{
                    cur_node.attr("stroke", null)
                    cur_node.attr("stroke-width", "0px")
                    removeValueFilter(cur_node.attr("fill"))
                }
           })
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
        <svg ref={svgRef} style={{cursor:"pointer"}}><g ref={gRef}></g></svg>
        {legend?.height == null &&  // only the "simple" legend has a hight attribute
        <Button variant="outlined" onClick={()=> {toggleUseVSUP()}}>Switch Encoding</Button>}
        <Button variant="outlined" onClick={()=> {clearFilter(colorSections)}}>Clear Filter</Button>
    </>
      
})
  
  