import { Box, Grid, Typography } from "@mui/material";
import { Button, Tooltip } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import React from "react";
import { AttributeSelectionTable } from "projection-space-explorer";
import { setPacoAttributes, setPacoConstraints } from "../../State/PacoSettingsDuck";
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import RotateLeftIcon from '@mui/icons-material/RotateLeft';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import SettingsIcon from '@mui/icons-material/Settings';
import { downloadImpl, map_shortname_to_smiles } from "../../Utility/Utils";
import * as d3v5 from "d3v5";

const mapStateToProps = (state: AppState) => ({
    pacoAttributes: state.pacoSettings?.pacoAttributes,
    pacoConstraints: state.pacoSettings?.pacoConstraints,
    pacoRef: state.pacoSettings?.pacoRef,
    dataset: state.dataset
});

const mapDispatchToProps = (dispatch) => ({
    setPacoAttributes: (atts) => dispatch(setPacoAttributes(atts)),
    setPacoConstraints: (consts) => dispatch(setPacoConstraints(consts)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  splitRef: any;
};

export const PacoTabPanel = connector(({setPacoAttributes, pacoAttributes, setPacoConstraints, pacoConstraints, pacoRef, dataset}: Props) => {
    let fileInput = React.useRef<any>();

    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <Typography variant="subtitle2" gutterBottom>
            Parallel Coordinates Settings
          </Typography>
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
            {/* TODO: also save chosen attributes? */}
            <AttributeSelectionTable attributes={pacoAttributes} setAttributes={setPacoAttributes} btnFullWidth={true}><SettingsIcon/>&nbsp;Choose Attributes</AttributeSelectionTable>
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
            <Tooltip title="Reset constraints to initial state">
                <Button fullWidth variant="outlined" aria-label="Reset constraints to initial state" onClick={() => {setPacoConstraints([...pacoConstraints])}}>
                    <RotateLeftIcon />&nbsp;Reset Constraints
                </Button>
            </Tooltip>
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
            <Grid container>
                <Grid item xs={6} paddingRight={1}>
                    <Tooltip title="Export constraints">
                        <Button fullWidth variant="outlined" color="primary" aria-label="Export constraints" onClick={() => {downloadConstraints(pacoRef.data[0].dimensions, dataset.columns)}}>
                            <FileDownloadIcon />&nbsp;Export
                        </Button>
                    </Tooltip>
                </Grid>
                <Grid item xs={6} paddingLeft={1}>
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
                        <Button fullWidth variant="outlined" color="primary" aria-label="Import constraints" onClick={() => fileInput.current.click()}>
                            <FileUploadIcon />&nbsp;Import
                        </Button>
                    </Tooltip>
                </Grid>
            </Grid>
        </Box>
      </div>
    );
  }
);

const downloadArrayAsCSV = (array, header) => {
    const csv_lines = array.map((row)=>{
        return Object.values(row).join(",")
    })
    let csv_content = header.join(",") + "\n"
    csv_content += csv_lines.join("\n")
    downloadImpl(csv_content, "parallel_coordinates_constraints.csv", "text/csv")
}

const downloadConstraints = (dimensions, columns) => {
    debugger;
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
                    let val = const_dimension.ticktext[n];
                    if(columns[const_dimension.label].metaInformation.imgSmiles){
                        val = map_shortname_to_smiles(val);
                    }
                    let constraint_object = { col: const_dimension.label, operator: "EQUALS", val1: val, val2: val }
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