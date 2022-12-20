import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import './PacoContext.scss';
import Plotly from 'plotly.js-dist';
import { Box, Button, IconButton, Tooltip, Typography } from '@mui/material';
import { useCancellablePromise } from 'projection-space-explorer';
import FileDownloadIcon from '@mui/icons-material/FileDownload';
import UpdateIcon from '@mui/icons-material/Update';
import FileUploadIcon from '@mui/icons-material/FileUpload';
import SettingsIcon from '@mui/icons-material/Settings';
// @ts-ignore
import { AttributeSelectionTable } from 'projection-space-explorer';
import { LoadingIndicatorView } from '../Overrides/Dataset/DatasetTabPanel';
import { ReactionCIMEBackendFromEnv } from '../Backend/ReactionCIMEBackend';
import { AppState } from '../State/Store';

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  triggerDatasetUpdate: state.handleDataset?.triggerUpdate,
});

const mapDispatchToProps = (dispatch) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};

const loading_area = 'loading_indicator_paco';
export const PacoContext = connector(function ({ dataset, triggerDatasetUpdate }: Props) {
  if (dataset == null || dataset.columns == null) return null;

  const paco_ref = React.useRef<any>();
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
    const constraint_dimensions = dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0);
    const all_constraints = [];
    for (const i in constraint_dimensions) {
      const const_dimension = constraint_dimensions[i];
      let constraintarray = const_dimension.constraintrange;
      if (!Array.isArray(constraintarray[0])) {
        // check, if it is a 1-dimensional array and transform it into a 2-d array
        constraintarray = [constraintarray];
      }
      for (const j in constraintarray) {
        const constraint = constraintarray[j];

        // handle numeric data
        if (dataset.columns[const_dimension.label].isNumeric) {
          const constraint_object = { col: const_dimension.label, operator: 'BETWEEN', val1: constraint[0], val2: constraint[1] };
          all_constraints.push(constraint_object);
        } else {
          // handle categorical data
          const lower = Math.ceil(constraint[0]);
          const upper = Math.floor(constraint[1]);
          for (let n = lower; n <= upper; n++) {
            // iterate over all real valued indices and add them to the constraints
            const constraint_object = { col: const_dimension.label, operator: 'EQUALS', val1: const_dimension.ticktext[n], val2: const_dimension.ticktext[n] };
            all_constraints.push(constraint_object);
          }
        }
      }
    }

    ReactionCIMEBackendFromEnv.updatePOIConstraints(dataset.info.path, all_constraints).then((res) => {
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
      const abort_controller = new AbortController(); // TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadPacoCSV function?
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

            function get_processed_info(values, key) {
              const current_constraints = constraints.filter((constraint) => constraint.col === key);

              // handle numeric data
              if (dataset.columns[key].isNumeric) {
                const val_range = dataset.columns[key].range.max - dataset.columns[key].range.min;
                const eps = val_range * 0.01;

                const constraintrange = [];
                for (const i in current_constraints) {
                  const constraint = current_constraints[i];
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
              const num_values = values.map((val) => distinct.indexOf(val));

              const constraintrange = [];
              for (const i in current_constraints) {
                const constraint = current_constraints[i];
                if (constraint.operator === 'EQUALS') {
                  const num_val = distinct.indexOf(constraint.val1);
                  constraintrange.push([num_val - 0.5, num_val + 0.5]); // have to add a small amount to get a range
                }
              }

              return { ticktext: distinct, values: num_values, tickvals: [...new Set(num_values)], constraintrange };
            }

            const cols = Object.keys(rows[0]);
            const dimensions = cols.map((v, i) => {
              const values = unpack(rows, v);
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

            Plotly.newPlot(paco_ref.current, [paco], layout, config);

            paco_ref.current
              .on('plotly_hover', function (data) {
                console.log('plotly_hover');
                console.log(data);
              })
              .on('plotly_unhover', function (data) {
                console.log('plotly_unhover');
                console.log(data);
              });
          },
          dataset.info.path,
          pacoShowColumns.map((col) => col.feature),
          cancellablePromise,
          abort_controller,
          loading_area,
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
                updateBackendConstraints(paco_ref.current.data[0].dimensions);
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
        <div style={{}} ref={paco_ref} />
      </Box>
      <LoadingIndicatorView area={loading_area} />
    </div>
  );
});
