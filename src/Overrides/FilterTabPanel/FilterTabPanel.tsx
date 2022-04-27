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
import { NOItemsInfo } from "../../Utility/NOItemsInfo";



const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  triggerDatasetUpdate: state.handleDataset?.triggerUpdate
});

const mapDispatchToProps = (dispatch) => ({
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
};

export const FilterTabPanel = connector(({dataset, triggerDatasetUpdate}: Props) => {

  // const [constraintsMap, setConstraintsMap] = React.useState({});
  const [constraints, setConstraints] = React.useState([]);
  const [constraintCols, setConstraintCols] = React.useState([]);
  let fileInput = React.useRef<any>();
  const [totalDataPoints, setTotalDataPoints] = React.useState(-1)

  
  React.useEffect(()=> {
    if(dataset != null){
      ReactionCIMEBackendFromEnv.loadPOIConstraints(dataset.info.path).then((res_constraints) => {
        const con_cols = [...new Set(res_constraints.map((con) => con.col))];

        setConstraintCols(con_cols)
        setConstraints(res_constraints)
      })
    }
  }, [dataset])

  return (dataset &&
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Typography variant="subtitle2" gutterBottom>
          Filter Settings
        </Typography>
        {/* TODO: always show info about NO items */}
        <NOItemsInfo variant="filterOutOfTotal"></NOItemsInfo> 
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
            }
          }}
        />}
      </Box>
      <Box paddingTop={1} paddingRight={2}>
          <FilterSettings 
            triggerDatasetUpdate={triggerDatasetUpdate}
            dataset={dataset}
            constraintCols={constraintCols}
            constraints={constraints}
            removeFilter={(col) => {
              setConstraintCols(constraintCols.filter((con) => con !== col))
            }}
          ></FilterSettings>
      </Box>

      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Grid container>
          <Grid item xs={6} paddingRight={1}>
            <Button
              fullWidth
              variant="outlined"
              onClick={() => {
                downloadConstraints(dataset.info.path)
              }}
            >
              <GetAppIcon />
              &nbsp;Export
            </Button>
          </Grid>

          <Grid item xs={6} paddingLeft={1}>
            <input
                style={{ display: 'none' }}
                accept={".csv"}
                ref={fileInput}
                // multiple
                type="file"
                onChange={(e) => {
                  uploadConstraints(e.target.files, dataset, triggerDatasetUpdate);
                }}
            />
            <Button
              fullWidth
              variant="outlined"
              onClick={() => {
                fileInput.current.click()
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


const downloadConstraints = (path) => {
  ReactionCIMEBackendFromEnv.downloadPOIConstraints(path);
}
const uploadConstraints = (files, dataset, triggerDatasetUpdate) => {
  if (files == null || files.length <= 0) {
      return;
  }
  var file = files[0];
  var fileName = file.name as string;

  if (fileName.endsWith("csv")) {
      ReactionCIMEBackendFromEnv.uploadPOIConstraints(dataset.info.path, file).then((res) => {
          if(res.msg === "ok" && triggerDatasetUpdate != null){
              triggerDatasetUpdate({
                  display: dataset.info.path,
                  path: dataset.info.path,
                  type: dataset.info.type,
                  uploaded: true
              })
          }
      })
  }
}