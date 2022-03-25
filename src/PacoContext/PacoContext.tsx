import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../State/Store";
import "./PacoContext.scss";
import Plotly from 'plotly.js-dist'
import * as d3v5 from "d3v5";
import { Button } from "@mui/material";
import { ReactionCIMEBackendFromEnv } from "../Backend/ReactionCIMEBackend";
import { useCancellablePromise } from "projection-space-explorer";
import { LoadingIndicatorView } from "../Overrides/Dataset/DatasetTabPanel";



 const mapStateToProps = (state: AppState) => ({
    dataset: state.dataset
  });
  
  const mapDispatchToProps = (dispatch) => ({
    // setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
  };
  
const loading_area = "loading_indicator_paco";
export const PacoContext = connector(function ({dataset}: Props) {
    if(dataset == null || dataset.columns == null)
        return null;

    let paco_ref = React.useRef<any>();
    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    // TODO: add user input to select, which columns to show
    const [pacoShowColumns, setPacoShowColumns] = React.useState(Object.keys(dataset.columns).filter((col) => dataset.columns[col].metaInformation.paco));

    const updateBackendConstraints = (dimensions) => {
        const constraint_dimensions = dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0)
        let all_constraints = []
        for(const i in constraint_dimensions){
            const const_dimension = constraint_dimensions[i]
            let constraintarray = const_dimension.constraintrange
            if(!Array.isArray(constraintarray[0])){ // check, if it is a 1-dimensional array and transform it into a 2-d array
                constraintarray = [constraintarray]
            }
            for(const j in constraintarray){
                const constraint = constraintarray[j]
                
                //handle numeric data
                if(dataset.columns[const_dimension.label].isNumeric){
                    let constraint_object = { col: const_dimension.label, operator: "BETWEEN", val1: constraint[0], val2: constraint[1] }
                    all_constraints.push(constraint_object)
                }else{ // handle categorical data
                    const lower = Math.ceil(constraint[0])
                    const upper = Math.floor(constraint[1])
                    for (var n = lower; n <= upper; n++) { // iterate over all real valued indices and add them to the constraints
                        let constraint_object = { col: const_dimension.label, operator: "EQUALS", val1: const_dimension.ticktext[n], val2: const_dimension.ticktext[n] }
                        all_constraints.push(constraint_object)
                    }
                }
                
            }
        }

        ReactionCIMEBackendFromEnv.updatePOIConstraints(dataset.info.path, all_constraints).then((res) => console.log(res))

    }

    React.useEffect(()=> {
        if(pacoShowColumns.length > 0){
            cancelPromises();
            let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadPacoCSV function?
            ReactionCIMEBackendFromEnv.loadPOIConstraints(dataset.info.path).then((constraints) => {
                ReactionCIMEBackendFromEnv.loadPacoCSV((rows) => {
                    function unpack(rows, key) {
                        return rows.map(function(row) {
                            let val = row[key];
                            if(dataset.columns[key].isNumeric)
                                val = parseFloat(val)
                            return val;
                        });
                    }
    
                    function get_processed_info(values, key){
                        const current_constraints = constraints.filter((constraint) => constraint.col === key);

                        //handle numeric data
                        if(dataset.columns[key].isNumeric){
                            const val_range = dataset.columns[key].range.max - dataset.columns[key].range.min
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

                        return {ticktext: distinct, values: num_values, tickvals: [...new Set(num_values)], constraintrange: constraintrange};
                    }
    
                    const cols = Object.keys(rows[0])
                    let dimensions = cols.map((v, i) => {
                        const values = unpack(rows, v)
                        const processed_info = get_processed_info(values, v);
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
    
                    // var color = unpack(rows, 'yield');
                    var paco = {
                        type: 'parcoords',
                        line: {
                            showscale: true,
                            reversescale: true,
                            // colorscale: 'YlOrRd',
                            // cmin: -4000,
                            // cmid: 0,
                            // cmax: -100,
                            // color: color,
                        },
                    
                        dimensions: dimensions
                    };
    
                    var layout = {
                        // width: 1500,
                        // height: 800, 
                        // hovermode: 'closest'
                    }
    
                    var config = {
                        responsive: true
                    }
    
                    Plotly.newPlot(paco_ref.current, [paco], layout, config);
    
                }, dataset.info.path, pacoShowColumns, cancellablePromise, abort_controller, loading_area)
            })
            
        }
        
        
    }, [dataset])
  
    return <div className="PacoParent">
        <Button onClick={() => {updateBackendConstraints(paco_ref.current.data[0].dimensions)}}>Update POIs</Button>
      <div style={{}} ref={paco_ref}></div>
      <LoadingIndicatorView 
        area={loading_area}
      />
    </div>
  })