import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import './PacoContext.scss';
import Plotly from 'plotly.js-dist';
import { Box } from '@mui/material';
import { selectVectors } from 'projection-space-explorer';
import { AppState } from '../State/Store';
import { arrayEquals, LIGHT_GREY, mapShortnameToSmiles, mapSmilesToShortname, RED } from '../Utility/Utils';
import { setPacoRef } from '../State/PacoSettingsDuck';

function unpack(col, rows, key) {
  return rows.map(function (row) {
    let val = row[key];
    if (col.isNumeric) val = parseFloat(val);
    return val;
  });
}

function getProcessedInfo(constraints, col, values, key) {
  const currentConstraints = constraints.filter((constraint) => constraint.col === key);

  // handle numeric data
  if (col.isNumeric) {
    const valRange = col.range.max - col.range.min;
    const eps = valRange * 0.01;

    const constraintrange = [];
    currentConstraints.forEach((constraint) => {
      if (constraint.operator === 'BETWEEN') {
        constraintrange.push([+constraint.val1, +constraint.val2]);
      } else if (constraint.operator === 'EQUALS') {
        constraintrange.push([+constraint.val1 - eps, +constraint.val1 + eps]); // have to add a small amount to get a range
      }
    });
    return { ticktext: undefined, values, tickvals: undefined, constraintrange };
  }

  // handle categorical data
  const distinct = [...new Set(values)];
  const numValues = values.map((val) => distinct.indexOf(val));

  const constraintrange = [];
  currentConstraints.forEach((constraint) => {
    if (constraint.operator === 'EQUALS') {
      const numVal = distinct.indexOf(constraint.val1);
      constraintrange.push([numVal - 0.5, numVal + 0.5]);
    }
  });

  if (col.metaInformation.imgSmiles) {
    return {
      ticktext: distinct.map((val) => mapSmilesToShortname(val as string)),
      values: numValues,
      tickvals: [...new Set(numValues)],
      constraintrange,
    };
  }

  return { ticktext: distinct, values: numValues, tickvals: [...new Set(numValues)], constraintrange };
}

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  pacoAttributes: state.pacoSettings?.pacoAttributes,
  pacoConstraints: state.pacoSettings?.pacoConstraints,
  currentAggregation: state.currentAggregation,
});

const mapDispatchToProps = (dispatch) => ({
  setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
  setPacoRef: (ref) => dispatch(setPacoRef(ref)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

// eslint-disable-next-line @typescript-eslint/no-shadow
export const PacoContext = connector(function ({ dataset, pacoAttributes, pacoConstraints, setPacoRef, setCurrentAggregation, currentAggregation }: Props) {
  if (dataset == null || dataset.columns == null) return null;

  const [pacoAggregation, setPacoAggregation] = React.useState([]);

  const line = {
    // colorscale: 'YlOrRd',
    // cmin: -4000,
    // cmid: 0,
    // cmax: -100,
    // color: color,
    showscale: false,
    reversescale: false,
    color: dataset.vectors.map((row) => 1),
    colorscale: [
      [0, LIGHT_GREY],
      [1, RED],
    ], // PSE_BLUE
  };
  const paco = {
    type: 'parcoords',
    line,
  };

  const layout = {
    padding: {
      top: 0,
      bottom: 0,
      left: 0,
      right: 0,
    },
    // width: 1500,
    // height: 800,
    // hovermode: 'closest'
  };

  const config = {
    responsive: true,
    displayModeBar: false,
  };

  const pacoRef = React.useRef<any>();

  React.useEffect(() => {
    setPacoRef(pacoRef?.current);
    // eslint-disable-next-line
  }, [pacoRef]);

  React.useEffect(() => {
    if (pacoAttributes != null) {
      const pacoShowColumns = pacoAttributes.filter((col) => col.show).map((value) => value.feature);
      if (pacoShowColumns.length > 0) {
        const cols = pacoShowColumns;
        const dimensions = cols.map((v, i) => {
          const values = unpack(dataset.columns[v], dataset.vectors, v);
          const processedInfo = getProcessedInfo(pacoConstraints, dataset.columns[v], values, v);
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

        // dimensions.push({ values: dataset.vectors.map((v) => v.__meta__.meshIndex), label: UNIQUE_ID, constraintrange: undefined, tickvals: [], ticktext: []})

        // var color = unpack(rows, 'yield');

        const newPaco = { ...paco, dimensions };
        Plotly.newPlot(pacoRef.current, [newPaco], layout, config);

        pacoRef.current.on('plotly_restyle', (data) => {
          // only change aggregation, if constraints were changed
          if (Object.keys(data[0]).filter((item) => item.includes('constraintrange')).length > 0) {
            // reset coloring of lines
            Plotly.restyle(pacoRef.current, { line: { ...line } }, [0]);

            const constraintDims = pacoRef.current.data[0].dimensions.filter((dim) => dim.constraintrange != null && dim.constraintrange.length > 0);
            if (constraintDims.length <= 0) {
              const agg = dataset.vectors.map((row) => row.__meta__.meshIndex);
              if (!arrayEquals(currentAggregation.aggregation, agg)) {
                setPacoAggregation(agg);
                setCurrentAggregation(agg);
              }
              return;
            }

            const filteredVectors = dataset.vectors.filter((row) => {
              let highlightItem = true;

              constraintDims.forEach((constDim) => {
                const col = constDim.label;
                const value = row[col];

                let constraintarray = constDim.constraintrange;
                if (!Array.isArray(constraintarray[0])) {
                  // check, if it is a 1-dimensional array and transform it into a 2-d array
                  constraintarray = [constraintarray];
                }

                let dimHighlightItem = false;
                constraintarray.forEach((constraint) => {
                  // handle numeric data
                  if (dataset.columns[col].isNumeric) {
                    dimHighlightItem = dimHighlightItem || (value > constraint[0] && value < constraint[1]);
                  } else {
                    // handle categorical data
                    const lower = Math.ceil(constraint[0]);
                    const upper = Math.floor(constraint[1]);
                    for (let n = lower; n <= upper; n++) {
                      // iterate over all real valued indices and add them to the constraints
                      let val = constDim.ticktext[n];
                      if (dataset.columns[constDim.label].metaInformation.imgSmiles) {
                        val = mapShortnameToSmiles(val);
                        dimHighlightItem = dimHighlightItem || value === val;
                      }
                    }
                  }
                });

                highlightItem = highlightItem && dimHighlightItem;
              });

              return highlightItem;
            });

            const agg = filteredVectors.map((row) => row.__meta__.meshIndex);
            if (!arrayEquals(currentAggregation.aggregation, agg)) {
              setPacoAggregation(agg);
              setCurrentAggregation(agg);
            }
          }
        });

        // const ticks = d3v5.select(paco_ref.current).selectAll(".tick text")
        // console.log(ticks)
        // ticks.on("mouseenter", (d) => {
        //     console.log("tick")
        //     console.log(d)
        // })
      }
    }
    // eslint-disable-next-line
  }, [dataset, pacoAttributes, pacoConstraints]);

  React.useEffect(() => {
    if (currentAggregation.aggregation != null && currentAggregation.aggregation.length > 0) {
      if (!arrayEquals(currentAggregation.aggregation, pacoAggregation)) {
        // const dims = paco_ref.current.data[0].dimensions;
        const color = dataset.vectors.map((row) => (currentAggregation.aggregation.includes(row.__meta__.meshIndex) ? 1 : 0));
        const newLine = { ...line, color };

        // var update = {
        //     line: new_line
        // }
        // Plotly.restyle(paco_ref.current, update, [0]);

        // use this to also reset constraints of paco
        const dimensions = pacoRef.current.data[0].dimensions.map((dim) => {
          return { ...dim, constraintrange: [] };
        });
        const newPaco = { ...paco, dimensions, line: newLine };
        Plotly.react(pacoRef.current, [newPaco], layout, config);

        setPacoAggregation([...currentAggregation.aggregation]);
      }
    } else {
      // reset coloring of lines
      Plotly.restyle(pacoRef.current, { line: { ...line } }, [0]);
    }
    // eslint-disable-next-line
  }, [currentAggregation]);

  return (
    <div className="PacoParent">
      {/* <div id="paco_tooltip" style={{position: "absolute", opacity: "0"}}></div> */}
      <Box style={{ clear: 'both' }} paddingLeft={2} paddingTop={0} paddingRight={0}>
        <div style={{}} ref={pacoRef} />
      </Box>
    </div>
  );
});
