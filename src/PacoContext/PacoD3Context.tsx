import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { Box, Button, IconButton, Tooltip } from '@mui/material';
import { AttributeSelectionTable, selectVectors } from 'projection-space-explorer';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import RotateLeftIcon from '@mui/icons-material/RotateLeft';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import SettingsIcon from '@mui/icons-material/Settings';
import * as d3v5 from 'd3v5';

import 'parcoord-es/dist/parcoords.css';
import ParCoords from 'parcoord-es';
import { downloadImpl } from '../Utility/Utils';
import { AppState } from '../State/Store';

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
  constraintDimensions.forEach((constDimension) => {
    let constraintarray = constDimension.constraintrange;
    if (!Array.isArray(constraintarray[0])) {
      // check, if it is a 1-dimensional array and transform it into a 2-d array
      constraintarray = [constraintarray];
    }
    constraintarray.forEach((constraint) => {
      // handle numeric data
      if (columns[constDimension.label].isNumeric) {
        const constraintObject = { col: constDimension.label, operator: 'BETWEEN', val1: constraint[0], val2: constraint[1] };
        allConstraints.push(constraintObject);
      } else {
        // handle categorical data
        const lower = Math.ceil(constraint[0]);
        const upper = Math.floor(constraint[1]);
        for (let n = lower; n <= upper; n++) {
          // iterate over all real valued indices and add them to the constraints
          const constraintObject = { col: constDimension.label, operator: 'EQUALS', val1: constDimension.ticktext[n], val2: constDimension.ticktext[n] };
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
  dataset: state.dataset,
});

const mapDispatchToProps = (dispatch) => ({
  setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};

export const PacoContext = connector(function ({ dataset }: Props) {
  if (dataset == null || dataset.columns == null || dataset.vectors == null) return null;

  const pacoRef = React.useRef<any>();
  const fileInput = React.useRef<any>();
  // TODO: also save chosen attributes?
  const [pacoAttributes, setPacoAttributes] = React.useState(
    Object.keys(dataset.columns).map((col) => {
      return { feature: col, show: dataset.columns[col].metaInformation.paco };
    }),
  );
  const [pacoConstraints, setPacoConstraints] = React.useState([]);

  React.useEffect(() => {
    const pacoShowColumns = pacoAttributes.filter((col) => col.show).map((value) => value.feature);
    if (pacoShowColumns.length > 0) {
      const cols = pacoShowColumns;
      //     const data = dataset.vectors.map((row) => {
      //         // TODO: filter out columns that should not be shown
      //         return row
      //     });
      const data = dataset.vectors;
      // const data = [{axis1: 5, axis2: 9, axis3: 4}, {axis1: 6, axis2: 1, axis3: 2}]
      const parcoords = ParCoords()(pacoRef.current);
      parcoords
        .mode('queue') // mode: "queue" --> for large dataset such that everything is still responsive during rendering
        .data(data)
        // .alpha(0.3)
        // .color(color)
        // .axisDots(0.2) // not working?
        // .hideAxis(["yield"])
        // .composite("darker")
        .render()
        .shadows()
        .reorderable()
        .brushMode('1D-axes'); // enable brushing
    }
  }, [dataset, pacoAttributes, pacoConstraints]);

  return (
    <div className="PacoParent">
      <div id="paco_tooltip" style={{ position: 'absolute', opacity: '0' }} />
      <Box style={{ clear: 'both' }} paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Box style={{ float: 'left' }} paddingTop={1}>
          <AttributeSelectionTable attributes={pacoAttributes} setAttributes={setPacoAttributes} btnFullWidth={false}>
            <SettingsIcon />
            &nbsp;Choose Attributes
          </AttributeSelectionTable>
        </Box>
        <Box style={{ float: 'right' }}>
          <Tooltip title="Reset constraints to initial state">
            <Button
              style={{ paddingRight: 25 }}
              variant="outlined"
              aria-label="Update Points of Interest"
              onClick={() => {
                setPacoConstraints([...pacoConstraints]);
              }}
            >
              <RotateLeftIcon />
              &nbsp;Reset Constraints
            </Button>
          </Tooltip>
          <Tooltip title="Export constraints">
            <IconButton
              color="primary"
              aria-label="Export POI constraints"
              onClick={() => {
                downloadConstraints(pacoRef.current.data[0].dimensions, dataset.columns);
              }}
            >
              <FileDownloadIcon />
            </IconButton>
          </Tooltip>
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
            <IconButton color="primary" aria-label="Import POI constraints" onClick={() => fileInput.current.click()}>
              <FileUploadIcon />
            </IconButton>
          </Tooltip>
        </Box>
      </Box>
      <Box style={{ clear: 'both' }} paddingLeft={2} paddingTop={1} paddingRight={2}>
        <div style={{}} ref={pacoRef} />
      </Box>
    </div>
  );
});
