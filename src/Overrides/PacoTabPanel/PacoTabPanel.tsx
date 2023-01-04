import { Box, Grid, Typography, Button, Tooltip } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import React from 'react';
import { AttributeSelectionTable } from 'projection-space-explorer';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import RotateLeftIcon from '@mui/icons-material/RotateLeft';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import SettingsIcon from '@mui/icons-material/Settings';
import * as d3v5 from 'd3v5';
import { downloadImpl, mapShortnameToSmiles } from '../../Utility/Utils';
import { PacoActions } from '../../State/PacoSettingsDuck';
import { AppState } from '../../State/Store';

const downloadArrayAsCSV = (array, header) => {
  const csvLines = array.map((row) => {
    return Object.values(row).join(',');
  });
  let csvContent = `${header.join(',')}\n`;
  csvContent += csvLines.join('\n');
  downloadImpl(csvContent, 'parallel_coordinates_constraints.csv', 'text/csv');
};

const downloadConstraints = (dimensions, columns) => {
  const constraintDimensions = dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0);
  if (constraintDimensions.length <= 0) return;

  const allConstraints = [];
  constraintDimensions.forEach((const_dimension) => {
    let constraintarray = const_dimension.constraintrange;
    if (!Array.isArray(constraintarray[0])) {
      // check, if it is a 1-dimensional array and transform it into a 2-d array
      constraintarray = [constraintarray];
    }
    constraintarray.forEach((constraint) => {
      // handle numeric data
      if (columns[const_dimension.label].isNumeric) {
        const constraintObject = { col: const_dimension.label, operator: 'BETWEEN', val1: constraint[0], val2: constraint[1] };
        allConstraints.push(constraintObject);
      } else {
        // handle categorical data
        const lower = Math.ceil(constraint[0]);
        const upper = Math.floor(constraint[1]);
        for (let n = lower; n <= upper; n++) {
          // iterate over all real valued indices and add them to the constraints
          let val = const_dimension.ticktext[n];
          if (columns[const_dimension.label].metaInformation.imgSmiles) {
            val = mapShortnameToSmiles(val);
          }
          const constraintObject = { col: const_dimension.label, operator: 'EQUALS', val1: val, val2: val };
          allConstraints.push(constraintObject);
        }
      }
    });
  });

  downloadArrayAsCSV(allConstraints, Object.keys(allConstraints[0]));
};
const uploadConstraints = (files, setConstraints) => {
  if (files == null || files.length <= 0) {
    return;
  }
  const file = files[0];

  const fileReader = new FileReader();
  fileReader.onload = (e) => {
    const data = d3v5.csvParse(e.target.result as string);
    setConstraints(data);
  };
  fileReader.readAsBinaryString(file);
};

const mapStateToProps = (state: AppState) => ({
  pacoAttributes: state.pacoSettings?.pacoAttributes,
  pacoConstraints: state.pacoSettings?.pacoConstraints,
  pacoRef: state.pacoSettings?.pacoRef,
  dataset: state.dataset,
});

const mapDispatchToProps = (dispatch) => ({
  setPacoAttributes: (atts) => dispatch(PacoActions.setPacoAttributes(atts)),
  setPacoConstraints: (consts) => dispatch(PacoActions.setPacoConstraints(consts)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  splitRef: any;
};

export const PacoTabPanel = connector(({ setPacoAttributes, pacoAttributes, setPacoConstraints, pacoConstraints, pacoRef, dataset }: Props) => {
  const fileInput = React.useRef<any>();

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Typography variant="subtitle2" gutterBottom>
          Parallel Coordinates Settings
        </Typography>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        {/* TODO: also save chosen attributes? */}
        <AttributeSelectionTable attributes={pacoAttributes} setAttributes={setPacoAttributes}></AttributeSelectionTable>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Tooltip title="Reset constraints to initial state">
          <Button
            fullWidth
            variant="outlined"
            aria-label="Reset constraints to initial state"
            onClick={() => {
              setPacoConstraints([...pacoConstraints]);
            }}
          >
            <RotateLeftIcon />
            &nbsp;Reset Constraints
          </Button>
        </Tooltip>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Grid container>
          <Grid item xs={6} paddingRight={1}>
            <Tooltip title="Export constraints">
              <Button
                fullWidth
                variant="outlined"
                color="primary"
                aria-label="Export constraints"
                onClick={() => {
                  downloadConstraints(pacoRef.data[0].dimensions, dataset.columns);
                }}
              >
                <FileDownloadIcon />
                &nbsp;Export
              </Button>
            </Tooltip>
          </Grid>
          <Grid item xs={6} paddingLeft={1}>
            <input
              style={{ display: 'none' }}
              accept=".csv"
              ref={fileInput}
              type="file"
              onChange={(e) => {
                uploadConstraints(e.target.files, setPacoConstraints);
              }}
            />
            <Tooltip title="Import constraints">
              <Button fullWidth variant="outlined" color="primary" aria-label="Import constraints" onClick={() => fileInput.current.click()}>
                <FileUploadIcon />
                &nbsp;Import
              </Button>
            </Tooltip>
          </Grid>
        </Grid>
      </Box>
    </div>
  );
});
