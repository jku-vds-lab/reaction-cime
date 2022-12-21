import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../State/Store";
import "./PacoContext.scss";
import Plotly from 'plotly.js-dist'
import { Box } from "@mui/material";
import { selectVectors } from "projection-space-explorer";
import { arrayEquals, LIGHT_GREY, map_shortname_to_smiles, map_smiles_to_shortname, RED } from "../Utility/Utils";
import { setPacoRef } from "../State/PacoSettingsDuck";

function unpack(col, rows, key) {
    return rows.map(function(row) {
        let val = row[key];
        if(col.isNumeric)
            val = parseFloat(val)
        return val;
    });
}

function get_processed_info(constraints, col, values, key){
    const current_constraints = constraints.filter((constraint) => constraint.col === key);

    //handle numeric data
    if(col.isNumeric){
        const val_range = col.range.max -col.range.min
        const eps = val_range*0.01

        let constraintrange = []
        for(const i in current_constraints){
            const constraint = current_constraints[i];
            if(constraint.operator === "BETWEEN"){
                constraintrange.push([+constraint.val1, +constraint.val2])
            }else if(constraint.operator === "EQUALS"){
                constraintrange.push([+constraint.val1-eps, +constraint.val1+eps]) // have to add a small amount to get a range
            }
        }
        return {ticktext: undefined, values: values, tickvals: undefined, constraintrange: constraintrange};
    }

    // handle categorical data
    const distinct = [...new Set(values)];
    const num_values = values.map((val) => distinct.indexOf(val))

    let constraintrange = []
    for(const i in current_constraints){
        const constraint = current_constraints[i];
        if(constraint.operator === "EQUALS"){
            const num_val = distinct.indexOf(constraint.val1)
            constraintrange.push([num_val-0.5, num_val+0.5]) // have to add a small amount to get a range
        }
    }

    if(col.metaInformation.imgSmiles){
        return {ticktext: distinct.map((val) => map_smiles_to_shortname(val as string)), values: num_values, tickvals: [...new Set(num_values)], constraintrange: constraintrange};
    }

    return {ticktext: distinct, values: num_values, tickvals: [...new Set(num_values)], constraintrange: constraintrange};
}



 const mapStateToProps = (state: AppState) => ({
    dataset: state.dataset,
    pacoAttributes: state.pacoSettings?.pacoAttributes,
    pacoConstraints: state.pacoSettings?.pacoConstraints,
    currentAggregation: state.currentAggregation,
  });
  
  const mapDispatchToProps = (dispatch) => ({
    setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
    setPacoRef: (ref) => dispatch(setPacoRef(ref)),
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
  };
  
export const PacoContext = connector(function ({dataset, pacoAttributes, pacoConstraints, setPacoRef, setCurrentAggregation, currentAggregation}: Props) {
    if(dataset == null || dataset.columns == null)
        return null;


    const [pacoAggregation, setPacoAggregation] = React.useState([]);

    const line = {
        // colorscale: 'YlOrRd',
        // cmin: -4000,
        // cmid: 0,
        // cmax: -100,
        // color: color,
        showscale: false,
        reversescale: false,
        color: dataset.vectors.map((row) => 1),
        colorscale: [[0, LIGHT_GREY], [1, RED]] //PSE_BLUE
    }
    const paco = {
        type: 'parcoords',
        line: line
    };
    
    const layout = {
        padding: {
            top: 0,
            bottom: 0,
            left: 0,
            right: 0,
        }
        // width: 1500,
        // height: 800, 
        // hovermode: 'closest'
    }

    const config = {
        responsive: true,
        displayModeBar: false
    }

    let paco_ref = React.useRef<any>();

    React.useEffect(() => {
        setPacoRef(paco_ref?.current)
    // eslint-disable-next-line
    }, [paco_ref])

    React.useEffect(()=> {
        if(pacoAttributes != null){
            const pacoShowColumns = pacoAttributes.filter((col) => col.show).map((value) => value.feature);
            if(pacoShowColumns.length > 0){
                const cols = pacoShowColumns;
                let dimensions = cols.map((v, i) => {
                    const values = unpack(dataset.columns[v], dataset.vectors, v)
                    const processed_info = get_processed_info(pacoConstraints, dataset.columns[v], values, v);
                    return {
                        values: processed_info.values,
                        label: v,
                        // multiselect: false,
                        constraintrange: processed_info.constraintrange,
                        // range: [Math.min(...values), Math.max(...values)],
                        tickvals: processed_info.tickvals,
                        ticktext: processed_info.ticktext,
                        // tickformat: ..., // https://plotly.com/javascript/reference/parcoords/#parcoords-dimensions-items-dimension-tickformat
                        // visible: true, //TODO: set to false, if, for example, datatype not recognized
                    }
                });

                // dimensions.push({ values: dataset.vectors.map((v) => v.__meta__.meshIndex), label: UNIQUE_ID, constraintrange: undefined, tickvals: [], ticktext: []})

                // var color = unpack(rows, 'yield');
                

                const new_paco = {...paco, dimensions: dimensions};
                Plotly.newPlot(paco_ref.current, [new_paco], layout, config);

                paco_ref.current.on("plotly_restyle", (data) => {
                    // only change aggregation, if constraints were changed
                    if(Object.keys(data[0]).filter((item) => item.includes("constraintrange")).length > 0){
                        
                        // reset coloring of lines
                        Plotly.restyle(paco_ref.current, {line: {...line}}, [0]);

                        const constraint_dims = paco_ref.current.data[0].dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0);
                        if(constraint_dims.length <= 0){
                            let agg = dataset.vectors.map((row) => row.__meta__.meshIndex);
                            if(!arrayEquals(currentAggregation.aggregation, agg)){
                                setPacoAggregation(agg)
                                setCurrentAggregation(agg);
                            }
                            return;
                        }

                        let filtered_vectors = dataset.vectors.filter((row) => {
                            let highlightItem = true;
                            for (const i in constraint_dims) {
                                const const_dim = constraint_dims[i];
                                const col = const_dim.label;
                                const value = row[col];

                                let constraintarray = const_dim.constraintrange
                                if(!Array.isArray(constraintarray[0])){ // check, if it is a 1-dimensional array and transform it into a 2-d array
                                    constraintarray = [constraintarray]
                                }

                                let dimHighlightItem = false;
                                for(const j in constraintarray){
                                    const constraint = constraintarray[j]
                                    
                                    //handle numeric data
                                    if(dataset.columns[col].isNumeric){
                                        dimHighlightItem = dimHighlightItem || (value > constraint[0] && value < constraint[1])
                                    }else{ // handle categorical data
                                        const lower = Math.ceil(constraint[0])
                                        const upper = Math.floor(constraint[1])
                                        for (var n = lower; n <= upper; n++) { // iterate over all real valued indices and add them to the constraints
                                            let val = const_dim.ticktext[n];
                                            if(dataset.columns[const_dim.label].metaInformation.imgSmiles){
                                                val = map_shortname_to_smiles(val);
                                                dimHighlightItem = dimHighlightItem || value === val
                                            }
                                        }
                                    }
                                    
                                }

                                highlightItem = highlightItem && dimHighlightItem;
                            }

                            return highlightItem;
                        });

                        const agg = filtered_vectors.map((row) => row.__meta__.meshIndex)
                        if(!arrayEquals(currentAggregation.aggregation, agg)){
                            setPacoAggregation(agg)
                            setCurrentAggregation(agg);
                        }

                    }
                    
                    
                });

                // const ticks = d3v5.select(paco_ref.current).selectAll(".tick text")
                // console.log(ticks)
                // ticks.on("mouseenter", (d) => {
                //     console.log("tick")
                //     console.log(d)
                // })
            }
        }
    // eslint-disable-next-line
    }, [dataset, pacoAttributes, pacoConstraints])

    React.useEffect(() => {
        if(currentAggregation.aggregation != null && currentAggregation.aggregation.length > 0){
            if(!arrayEquals(currentAggregation.aggregation, pacoAggregation)){

                // const dims = paco_ref.current.data[0].dimensions;
                const color = dataset.vectors.map((row) => currentAggregation.aggregation.includes(row.__meta__.meshIndex) ? 1 : 0);
                const new_line = {...line, color: color};
                
                // var update = {
                //     line: new_line
                // }
                // Plotly.restyle(paco_ref.current, update, [0]);

                // use this to also reset constraints of paco
                const dimensions = paco_ref.current.data[0].dimensions.map((dim) => { return {...dim, constraintrange: []}});
                const new_paco = {...paco, dimensions: dimensions, line: new_line};
                Plotly.react(paco_ref.current, [new_paco], layout, config);

                setPacoAggregation([...currentAggregation.aggregation])

                return;
            }
        }
        else{
            // reset coloring of lines
            Plotly.restyle(paco_ref.current, {line: {...line}}, [0]);
        }
    // eslint-disable-next-line
    }, [currentAggregation])
  
    return <div className="PacoParent">
        {/* <div id="paco_tooltip" style={{position: "absolute", opacity: "0"}}></div> */}
        <Box style={{clear: "both"}} paddingLeft={2} paddingTop={0} paddingRight={0}>
            <div style={{}} ref={paco_ref}></div>
        </Box>
    </div>
  })