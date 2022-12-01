import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../State/Store";
import { Box, Button, IconButton, Tooltip } from "@mui/material";
import { AttributeSelectionTable, selectVectors } from "projection-space-explorer";
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import RotateLeftIcon from '@mui/icons-material/RotateLeft';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import SettingsIcon from '@mui/icons-material/Settings';
import { downloadImpl } from "../Utility/Utils";
import * as d3v5 from "d3v5";

import 'parcoord-es/dist/parcoords.css';
import ParCoords from 'parcoord-es';

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

    return {ticktext: distinct, values: num_values, tickvals: [...new Set(num_values)], constraintrange: constraintrange};
}

const downloadArrayAsCSV = (array, header) => {
    const csv_lines = array.map((row)=>{
        return Object.values(row).join(",")
    })
    let csv_content = header.join(",") + "\n"
    csv_content += csv_lines.join("\n")
    downloadImpl(csv_content, "parallel_coordinates_constraints.csv", "text/csv")
}

const downloadConstraints = (dimensions, columns) => {
    const constraint_dimensions = dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0)
    if(constraint_dimensions.length <= 0)
        return;

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
            if(columns[const_dimension.label].isNumeric){
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
    downloadArrayAsCSV(all_constraints, Object.keys(all_constraints[0]))
}
const uploadConstraints = (files, setConstraints) => {
    if (files == null || files.length <= 0) {
        return;
    }
    var file = files[0];

    const fileReader = new FileReader()
    fileReader.onload = (e) => {
        const data = d3v5.csvParse(e.target.result as string)
        setConstraints(data)
    }
    fileReader.readAsBinaryString(file)
}

 const mapStateToProps = (state: AppState) => ({
    dataset: state.dataset,
  });
  
  const mapDispatchToProps = (dispatch) => ({
    setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
  };
  
export const PacoContext = connector(function ({dataset}: Props) {
    if(dataset == null || dataset.columns == null || dataset.vectors == null)
        return null;

    let paco_ref = React.useRef<any>();
    let fileInput = React.useRef<any>();
    // TODO: also save chosen attributes?
    const [pacoAttributes, setPacoAttributes] = React.useState(Object.keys(dataset.columns).map((col) => {
        return {feature: col, show: dataset.columns[col].metaInformation.paco};
    }));
    const [pacoConstraints, setPacoConstraints] = React.useState([])

    React.useEffect(()=> {
        const pacoShowColumns = pacoAttributes.filter((col) => col.show).map((value) => value.feature);
        if(pacoShowColumns.length > 0){
            const cols = pacoShowColumns;
        //     const data = dataset.vectors.map((row) => {
        //         // TODO: filter out columns that should not be shown
        //         return row
        //     });
            const data = dataset.vectors;
            // const data = [{axis1: 5, axis2: 9, axis3: 4}, {axis1: 6, axis2: 1, axis3: 2}]
            var parcoords = ParCoords()(paco_ref.current);
            parcoords
                .mode("queue") // mode: "queue" --> for large dataset such that everything is still responsive during rendering
                .data(data)
                // .alpha(0.3)
                // .color(color)
                // .axisDots(0.2) // not working?
                // .hideAxis(["yield"])
                // .composite("darker")
                .render()
                .shadows()
                .reorderable()
                .brushMode("1D-axes");  // enable brushing
        }
    }, [dataset, pacoAttributes, pacoConstraints])
  
    return <div className="PacoParent">
        <div id="paco_tooltip" style={{position: "absolute", opacity: "0"}}></div>
        <Box style={{clear: "both"}} paddingLeft={2} paddingTop={1} paddingRight={2}>
            <Box style={{float: "left"}} paddingTop={1}>
                <AttributeSelectionTable attributes={pacoAttributes} setAttributes={setPacoAttributes} btnFullWidth={false}><SettingsIcon/>&nbsp;Choose Attributes</AttributeSelectionTable>
            </Box>
            <Box style={{float: "right"}}>
                <Tooltip title="Reset constraints to initial state">
                    <Button style={{paddingRight: 25}} variant="outlined" aria-label="Update Points of Interest" onClick={() => {setPacoConstraints([...pacoConstraints])}}>
                        <RotateLeftIcon />&nbsp;Reset Constraints
                    </Button>
                </Tooltip>
                <Tooltip title="Export constraints">
                    <IconButton color="primary" aria-label="Export POI constraints" onClick={() => {downloadConstraints(paco_ref.current.data[0].dimensions, dataset.columns)}}>
                        <FileDownloadIcon />
                    </IconButton>
                </Tooltip>
                <input
                        style={{ display: 'none' }}
                        accept={".csv"}
                        ref={fileInput}
                        type="file"
                        onChange={(e) => {
                            uploadConstraints(e.target.files, setPacoConstraints);
                        }}
                    />
                <Tooltip title="Import constraints">
                    <IconButton color="primary" aria-label="Import POI constraints" onClick={() => fileInput.current.click()}>
                        <FileUploadIcon />
                    </IconButton>
                </Tooltip>
            </Box>
        </Box>
        <Box style={{clear: "both"}} paddingLeft={2} paddingTop={1} paddingRight={2}>
            <div style={{}} ref={paco_ref}></div>
        </Box>
    </div>
  })