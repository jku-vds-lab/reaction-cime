import { Grid, Typography } from "@mui/material";
import { Box, Button } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import GetAppIcon from "@mui/icons-material/GetApp";
import { AppState } from "../../State/Store";
import React from "react";
// @ts-ignore
import { SelectFeatureComponent } from "projection-space-explorer";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { FilterSettings } from "./FilterSettings";
import FileUploadIcon from '@mui/icons-material/FileUpload';
import FilterAltIcon from '@mui/icons-material/FilterAlt';



const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
});

const mapDispatchToProps = (dispatch) => ({
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
};

export const FilterTabPanel = connector(({dataset}: Props) => {

  // const [constraintsMap, setConstraintsMap] = React.useState({});
  const [constraints, setConstraints] = React.useState([]);
  const [constraintCols, setConstraintCols] = React.useState([]);

  
  React.useEffect(()=> {
    if(dataset != null){
      ReactionCIMEBackendFromEnv.loadPOIConstraints(dataset.info.path).then((res_constraints) => {
        const con_cols = [...new Set(res_constraints.map((con) => con.col))];

        setConstraintCols(con_cols)
        setConstraints(res_constraints)
        // let tempConstraintsMap = {}
        // for (const i in constraints) {
        //   const constraint = constraints[i];
        //   if(Object.keys(tempConstraintsMap).includes(constraint.col)){
        //     tempConstraintsMap[constraint.col].push({operator: constraint.operator, val1: constraint.val1, val2: constraint.val2})
        //   }else{
        //     tempConstraintsMap[constraint.col] = [{operator: constraint.operator, val1: constraint.val1, val2: constraint.val2}]
        //   }
        // }
        // setConstraintsMap(()=> tempConstraintsMap)
      })
    }
  }, [dataset])

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Typography variant="subtitle2" gutterBottom>
          Filter Settings
        </Typography>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        {dataset && <SelectFeatureComponent
          column_info={dataset.columns}
          label="filter"
          default_val={undefined}
          // categoryOptions={Object.keys(dataset.columns).filter((col) => !Object.keys(constraintsMap).includes(col))}
          categoryOptions={Object.keys(dataset.columns).filter((col) => !constraintCols.includes(col))}
          onChange={(newValue) => {
            if(Object.keys(dataset.columns).includes(newValue)){
              setConstraintCols([...constraintCols, newValue]);
              // if(dataset.columns[newValue].isNumeric){
              //   let tempConstraintsMap = {...constraintsMap}
              //   tempConstraintsMap[newValue] = [{operator: "BETWEEN", val1: undefined, val2: undefined}] // defaults to min and max
              //   setConstraintsMap(tempConstraintsMap)
              // }else{
              //     let tempConstraintsMap = {...constraintsMap}
              //     tempConstraintsMap[newValue] = []
              //     setConstraintsMap(tempConstraintsMap)
              // }
            }
          }}
        />}
      </Box>
      <Box paddingTop={1} paddingRight={2}>
          <FilterSettings 
            dataset={dataset}
            constraintCols={constraintCols}
            constraints={constraints}
            // constraintsMap={constraintsMap}
            removeFilter={(col) => {
              
              setConstraintCols(constraintCols.filter((con) => con !== col))
              // let tempConstraintsMap = {...constraintsMap}
              // delete tempConstraintsMap[col];
              // setConstraintsMap(tempConstraintsMap)
            }}
          ></FilterSettings>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Button
          fullWidth
          variant="outlined"
          onClick={() => {
          }}
        >
          <FilterAltIcon />
          &nbsp;Apply Filter
        </Button>
      </Box>

      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Grid container>
          <Grid item xs={6} paddingRight={1}>
            <Button
              fullWidth
              variant="outlined"
              onClick={() => {
                console.log("TODO: export filter")
              }}
            >
              <GetAppIcon />
              &nbsp;Export
            </Button>
          </Grid>

          <Grid item xs={6} paddingLeft={1}>
            <Button
              fullWidth
              variant="outlined"
              onClick={() => {
                console.log("TODO: import filter")
              }}
            >
              <FileUploadIcon />
              &nbsp;Import
            </Button>
          </Grid>
        </Grid>
      </Box>
    </div>
  );
});
