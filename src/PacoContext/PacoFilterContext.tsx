import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import './PacoContext.scss';
import Plotly from 'plotly.js-dist';
import { Box, Button, IconButton, Tooltip, Typography } from '@mui/material';
import { useCancellablePromise } from 'projection-space-explorer';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import UpdateIcon from '@mui/icons-material/Update';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import { LoadingIndicatorView } from '../Overrides/Dataset/DatasetTabPanel';
import { ReactionCIMEBackendFromEnv } from '../Backend/ReactionCIMEBackend';
import { AppState } from '../State/Store';

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  triggerDatasetUpdate: state.handleDataset?.triggerUpdate,
});

const connector = connect(mapStateToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

const loadingArea = 'loading_indicator_paco';
export const PacoContext = connector(({ dataset, triggerDatasetUpdate }: Props) => {
  if (dataset == null || dataset.columns == null) return null;

  const pacoRef = React.useRef<any>();
  const fileInput = React.useRef<any>();
  const { cancellablePromise, cancelPromises } = useCancellablePromise();

  // const [pacoShowColumns, setPacoShowColumns] = React.useState(Object.keys(dataset.columns).filter((col) => dataset.columns[col].metaInformation.paco));
  const [pacoAttributes, setPacoAttributes] = React.useState(
    Object.keys(dataset.columns).map((col) => {
      return { feature: col, show: dataset.columns[col].metaInformation.paco };
    }),
  );
  const [totalDataPoints, setTotalDataPoints] = React.useState(-1);

  const updateBackendConstraints = (dimensions) => {
    const constraintDimensions = dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0);
    const allConstraints = [];
    for (let i = 0; i < constraintDimensions.length; i++) {
      const constDimension = constraintDimensions[i];
      let constraintarray = constDimension.constraintrange;
      if (!Array.isArray(constraintarray[0])) {
        // check, if it is a 1-dimensional array and transform it into a 2-d array
        constraintarray = [constraintarray];
      }
      for (let j = 0; j < constraintarray; j++) {
        const constraint = constraintarray[j];

        // handle numeric data
        if (dataset.columns[constDimension.label].isNumeric) {
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
      }
    }

    ReactionCIMEBackendFromEnv.updatePOIConstraints(dataset.info.path, allConstraints).then((res) => {
      if (res.msg === 'ok' && triggerDatasetUpdate != null) {
        triggerDatasetUpdate({
          display: dataset.info.path,
          path: dataset.info.path,
          type: dataset.info.type,
          uploaded: true,
        });
      }
    });
  };

  const downloadConstraints = () => {
    ReactionCIMEBackendFromEnv.downloadPOIConstraints(dataset.info.path);
  };
  const uploadConstraints = (files) => {
    if (files == null || files.length <= 0) {
      return;
    }
    const file = files[0];
    const fileName = file.name as string;

    if (fileName.endsWith('csv')) {
      ReactionCIMEBackendFromEnv.uploadPOIConstraints(dataset.info.path, file).then((res) => {
        if (res.msg === 'ok' && triggerDatasetUpdate != null) {
          triggerDatasetUpdate({
            display: dataset.info.path,
            path: dataset.info.path,
            type: dataset.info.type,
            uploaded: true,
          });
        }
      });
    }
  };

  React.useEffect(() => {
    const pacoShowColumns = pacoAttributes.filter((col) => col.show);
    if (pacoShowColumns.length > 0) {
      cancelPromises();
      const abortController = new AbortController(); // TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadPacoCSV function?
      ReactionCIMEBackendFromEnv.loadPOIConstraints(dataset.info.path).then((constraints) => {
        ReactionCIMEBackendFromEnv.loadPacoCSV(
          (rows) => {
            setTotalDataPoints(rows.length);

            function unpack(rows, key) {
              return rows.map(function (row) {
                let val = row[key];
                if (dataset.columns[key].isNumeric) val = parseFloat(val);
                return val;
              });
            }

            function getProcessedInfo(values, key) {
              const currentConstraints = constraints.filter((constraint) => constraint.col === key);

              // handle numeric data
              if (dataset.columns[key].isNumeric) {
                const valRange = dataset.columns[key].range.max - dataset.columns[key].range.min;
                const eps = valRange * 0.01;

                const constraintrange = [];
                for (let i = 0; i < currentConstraints; i++) {
                  const constraint = currentConstraints[i];
                  if (constraint.operator === 'BETWEEN') {
                    constraintrange.push([+constraint.val1, +constraint.val2]);
                  } else if (constraint.operator === 'EQUALS') {
                    constraintrange.push([+constraint.val1 - eps, +constraint.val1 + eps]); // have to add a small amount to get a range
                  }
                }
                return { ticktext: undefined, values, tickvals: undefined, constraintrange };
              }

              // handle categorical data
              const distinct = [...new Set(values)];
              const numValues = values.map((val) => distinct.indexOf(val));

              const constraintrange = [];
              for (let i = 0; i < currentConstraints; i++) {
                const constraint = currentConstraints[i];
                if (constraint.operator === 'EQUALS') {
                  const numVal = distinct.indexOf(constraint.val1);
                  constraintrange.push([numVal - 0.5, numVal + 0.5]); // have to add a small amount to get a range
                }
              }

              return { ticktext: distinct, values: numValues, tickvals: [...new Set(numValues)], constraintrange };
            }

            const cols = Object.keys(rows[0]);
            const dimensions = cols.map((v, i) => {
              const values = unpack(rows, v);
              const processedInfo = getProcessedInfo(values, v);
              return {
                values: processedInfo.values,
                label: v,
                // multiselect: false,
                constraintrange: processedInfo.constraintrange,
                // range: [Math.min(...values), Math.max(...values)],
                tickvals: processedInfo.tickvals,
                ticktext: processedInfo.ticktext,
                // tickformat: ..., // https://plotly.com/javascript/reference/parcoords/#parcoords-dimensions-items-dimension-tickformat
                // visible: true, //TODO: set to false, if, for example, datatype not recognized
              };
            });

            // var color = unpack(rows, 'yield');
            const paco = {
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

              dimensions,
            };

            const layout = {
              padding: {
                top: 0,
              },
              // width: 1500,
              // height: 800,
              // hovermode: 'closest'
            };

            const config = {
              responsive: true,
            };

            Plotly.newPlot(pacoRef.current, [paco], layout, config);

            pacoRef.current
              .on('plotly_hover', (data) => {
                console.log('plotly_hover');
                console.log(data);
              })
              .on('plotly_unhover', (data) => {
                console.log('plotly_unhover');
                console.log(data);
              });
          },
          dataset.info.path,
          pacoShowColumns.map((col) => col.feature),
          cancellablePromise,
          abortController,
          loadingArea,
        );
      });
    }
  }, [dataset, pacoAttributes]);

  return (
    <div className="PacoParent">
      <div id="paco_tooltip" style={{ position: 'absolute', opacity: '0' }} />
      <Box style={{ clear: 'both' }} paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Box style={{ float: 'left' }} paddingTop={1}>
          {/* <AttributeSelectionTable attributes={pacoAttributes} setAttributes={setPacoAttributes} btnFullWidth={false}><SettingsIcon/>&nbsp;Choose Attributes</AttributeSelectionTable> */}
        </Box>
        <Box style={{ float: 'right' }}>
          <Typography color="textSecondary" variant="body2">
            Currently showing {dataset.vectors.length} out of {totalDataPoints} items
          </Typography>
          <Tooltip title="Update Points of Interest">
            <Button
              style={{ paddingRight: 25 }}
              variant="outlined"
              aria-label="Update Points of Interest"
              onClick={() => {
                updateBackendConstraints(pacoRef.current.data[0].dimensions);
              }}
            >
              <UpdateIcon />
              &nbsp;Update POIs
            </Button>
          </Tooltip>
          <Tooltip title="Export POI constraints">
            <IconButton
              color="primary"
              aria-label="Export POI constraints"
              onClick={() => {
                downloadConstraints();
              }}
            >
              <FileDownloadIcon />
            </IconButton>
          </Tooltip>
          <input
            style={{ display: 'none' }}
            accept=".csv"
            ref={fileInput}
            // multiple
            type="file"
            onChange={(e) => {
              uploadConstraints(e.target.files);
            }}
          />
          <Tooltip title="Import POI constraints">
            <IconButton color="primary" aria-label="Import POI constraints" onClick={() => fileInput.current.click()}>
              <FileUploadIcon />
            </IconButton>
          </Tooltip>
        </Box>
      </Box>
      <Box style={{ clear: 'both' }} paddingLeft={2} paddingTop={1} paddingRight={2}>
        <div style={{}} ref={pacoRef} />
      </Box>
      <LoadingIndicatorView area={loadingArea} />
    </div>
  );
});
