import { Grid, Typography, Box, Button } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import GetAppIcon from '@mui/icons-material/GetApp';
import React from 'react';
import { toSentenceCase, SelectFeatureComponent } from 'projection-space-explorer';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import { AppState } from '../../State/Store';
// @ts-ignore
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { FilterSettings } from './FilterSettings';
import { NOItemsInfo } from '../../Utility/NOItemsInfo';
import { ExceptionSettings } from './ExceptionSettings';

const downloadConstraints = (path) => {
  ReactionCIMEBackendFromEnv.downloadPOIConstraints(path);
};
const uploadConstraints = (files, dataset, triggerDatasetUpdate, state) => {
  if (files == null || files.length <= 0) {
    return;
  }
  const file = files[0];
  const fileName = file.name as string;

  if (fileName.endsWith('csv')) {
    ReactionCIMEBackendFromEnv.uploadPOIConstraints(dataset.info.path, file).then((res) => {
      if (res.msg === 'ok' && triggerDatasetUpdate != null) {
        triggerDatasetUpdate(
          {
            display: dataset.info.path,
            path: dataset.info.path,
            type: dataset.info.type,
            uploaded: true,
          },
          state,
        );
      }
    });
  }
};

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  triggerDatasetUpdate: state.handleDataset?.triggerUpdate,
  state,
});

const mapDispatchToProps = (dispatch) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

export const FilterTabPanel = connector(({ dataset, triggerDatasetUpdate, state }: Props) => {
  // const [constraintsMap, setConstraintsMap] = React.useState({});
  const [constraints, setConstraints] = React.useState([]);
  const [constraintCols, setConstraintCols] = React.useState([]);
  const fileInput = React.useRef<any>();

  React.useEffect(() => {
    if (dataset != null) {
      ReactionCIMEBackendFromEnv.loadPOIConstraints(dataset.info.path).then((res_constraints) => {
        const conCols = [...new Set(res_constraints.map((con) => con.col))];
        setConstraintCols(conCols);
        setConstraints(res_constraints);
      });
    }
  }, [dataset]);

  return (
    dataset && (
      <div style={{ display: 'flex', flexDirection: 'column', height: '100%', overflowY: 'auto' }}>
        <Box paddingX={2} paddingTop={2} paddingBottom={1}>
          <Typography variant="subtitle2" gutterBottom>
            Filter settings
          </Typography>
          <Typography variant="body2" color="textSecondary" gutterBottom>
            Adjust filter settings to show a different subset of {state.globalLabels.itemLabelPlural} in the front end.{' '}
            {/* TODO: always show info about NO items */}
            <NOItemsInfo variant="filterOutOfTotal" />
          </Typography>
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <SelectFeatureComponent
            column_info={dataset.columns}
            label="filter"
            default_val={undefined}
            // categoryOptions={Object.keys(dataset.columns).filter((col) => !Object.keys(constraintsMap).includes(col))}
            categoryOptions={Object.keys(dataset.columns).filter((col) => !constraintCols.includes(col))}
            onChange={(newValue) => {
              if (Object.keys(dataset.columns).includes(newValue)) {
                setConstraintCols([...constraintCols, newValue]);
              }
            }}
          />
        </Box>
        <Box paddingTop={1} paddingRight={2}>
          <FilterSettings
            state={state}
            triggerDatasetUpdate={triggerDatasetUpdate}
            dataset={dataset}
            constraintCols={constraintCols}
            constraints={constraints}
            removeFilter={(col) => {
              setConstraintCols(constraintCols.filter((con) => con !== col));
            }}
          />
        </Box>
        {/* TODO: this downloads the current filter file. this file can later be uploaded again -> not really necessary because filters are saved in backend */}
        {/* <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <Grid container>
            <Grid item xs={6} paddingRight={1}>
              <Button
                fullWidth
                variant="outlined"
                onClick={() => {
                  downloadConstraints(dataset.info.path);
                }}
              >
                <GetAppIcon />
                &nbsp;Export
              </Button>
            </Grid>

            <Grid item xs={6} paddingLeft={1}>
              <input
                style={{ display: 'none' }}
                accept=".csv"
                ref={fileInput}
                // multiple
                type="file"
                onChange={(e) => {
                  uploadConstraints(e.target.files, dataset, triggerDatasetUpdate, state);
                }}
              />
              <Button
                fullWidth
                variant="outlined"
                onClick={() => {
                  fileInput.current.click();
                }}
              >
                <FileUploadIcon />
                &nbsp;Import
              </Button>
            </Grid>
          </Grid>
        </Box> */}
        <ExceptionSettings state={state} triggerDatasetUpdate={triggerDatasetUpdate} dataset={dataset} />
      </div>
    )
  );
});
