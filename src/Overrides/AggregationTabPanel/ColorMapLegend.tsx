import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { Button, InputLabel, MenuItem, Select } from "@mui/material";
import { D3_CONTINUOUS_COLOR_SCALE_LIST, setAggregateColorScale, addValueFilter, removeValueFilter, toggleUseVSUP, clearValueFilter, setAggregateColorMapScale } from "../../State/AggregateSettingsDuck";
import * as vsup from "vsup";
import * as d3 from 'd3v5';

const mapStateToProps = (state: AppState) => ({
    colorScale: state.aggregateSettings?.colorscale,
    useVSUP: state.aggregateSettings?.useVSUP,
    aggregateColor: state.aggregateSettings?.aggregateColor,
    valueFilter: state.aggregateSettings?.valueFilter,
    valueRange: state.aggregateSettings?.valueRange,
    uncertaintyRange: state.aggregateSettings?.uncertaintyRange,
});

const mapDispatchToProps = (dispatch) => ({
    setAggregateColorScale: value => dispatch(setAggregateColorScale(value)),
    toggleUseVSUP: () => dispatch(toggleUseVSUP()),
    addValueFilter: (value) => dispatch(addValueFilter(value)),
    removeValueFilter: (value) => dispatch(removeValueFilter(value)),
    clearValueFilter: () => dispatch(clearValueFilter()),
    setAggregateColorMapScale: (scale) => dispatch(setAggregateColorMapScale(scale)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    selectAttribute: {key:string, name:string, col_info:any}
};



export const ColorMapLegend = connector(({colorScale, setAggregateColorScale, selectAttribute, toggleUseVSUP, addValueFilter, removeValueFilter, clearValueFilter, aggregateColor, useVSUP, valueFilter, valueRange, uncertaintyRange, setAggregateColorMapScale }: Props) => {
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
    const [legend, setLegend] = React.useState(null);

    const clearFilter = (color_sections) => {
        clearValueFilter()
        
        color_sections.attr("stroke", null)
        color_sections.attr("stroke-width", "0px")
    }


    React.useEffect(() => {
        if(valueRange != null){
            // visit https://github.com/uwdata/vsup for more info about VSUP vs bivariate colorscale encoding
            var vDom = [valueRange.min, valueRange.max];
            let colLegend, scale, quantization;

            // simple encoding
            if(uncertaintyRange == null){
                scale = d3.scaleQuantize()
                    .domain(vDom)
                    .range(d3.quantize(d3[colorScale], 8));
                vDom = scale.domain()
                colLegend = vsup.legend.simpleLegend(scale).title(aggregateColor.value_col)
            }else{ // bivariate encoding
                const uDom = [uncertaintyRange.min, uncertaintyRange.max];
                if(useVSUP){
                    quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
                    scale = vsup.scale().quantize(quantization).range(d3[colorScale]);
                    colLegend = vsup.legend.arcmapLegend(scale);
                }else{
                    quantization = vsup.squareQuantization(4).valueDomain(vDom).uncertaintyDomain(uDom);
                    scale = vsup.scale().quantize(quantization).range(d3[colorScale]);
                    colLegend = vsup.legend.heatmapLegend(scale);
                }
                colLegend
                    .vtitle(aggregateColor.value_col)
                    .utitle(aggregateColor.uncertainty_col)
                    
            }
            
            // TODO: for extremely small values, we could add some scaling function that scales the values by some e-5 and add "e-5" to the label --> the ticks would then for example be "0.01" and the title "columnname e-5"
            // automatically obtain a suitable precision value
            const p = d3.precisionFixed((vDom[1]-vDom[0])/8);
            if(p > 2){
                colLegend.format(".2")
            }else if(p > 0){
                colLegend.format("." + p + "r");
            }
            
            setLegend(() => colLegend)
            setAggregateColorMapScale(scale)
        }
        // eslint-disable-next-line
    }, [useVSUP, valueRange, uncertaintyRange, colorScale, aggregateColor])

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
            .on("mouseover", (d, i) => {
                // TODO: add hover interaction --> highlight areas in background somehow (e.g. have a component that draws a border around the selected areas)
                legend_container.select(".hover_clone").remove();
                const hover_el = d3.select(color_sections.nodes()[i]).clone();
                hover_el.attr("class", "hover_clone");
                hover_el.attr("fill", "#1f77b4");
                
                // remove temporary hover element when we leave it
                hover_el.on("mouseout", () => {
                    hover_el.remove()
                })

                // when clicking, we want to select the underlying section
                hover_el.on("click", () => {
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
            })

            svgElement.on("mouseout", () => {
                legend_container.select(".hover_clone").remove();
            });
            
        }
    // eslint-disable-next-line
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
            <Button variant="outlined" onClick={()=> {toggleUseVSUP()}}>Switch Encoding</Button>
        }
        {(valueFilter != null && valueFilter.length > 0) && <Button variant="outlined" onClick={()=> {clearFilter(colorSections)}}>Clear Filter</Button>}
        
    </>
      
})
  
  