import * as React from 'react';
// adapted from PSEs CoralLegend
import { connect, ConnectedProps } from 'react-redux';
import { Handler } from 'vega-tooltip';
import { makeStyles } from '@mui/styles';
import './FeatureLegend.scss';
import { Table, TableBody, TableCell, TableHead, TableRow, Tooltip, Typography } from '@mui/material';
import * as vegaImport from 'vega';
import { DefaultLegend, FeatureType, IVector, RootState } from 'projection-space-explorer';
import { InfoOutlined } from '@mui/icons-material';
import VegaDensity from './VegaHelpers/VegaDensity';
import BarChart from './VegaHelpers/BarChart';
import VegaDate from './VegaHelpers/VegaDate';
import { mapSmilesToShortname } from '../../Utility/Utils';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';

export function formatSMILESTooltip(value: any, valueToHtml: (value: any) => string, maxDepth: number): string {
  let content = '';

  content += '<table>';
  content += `<tr><td class="key">ratio:</td><td class="value">${valueToHtml(value.feature)}</td></tr>`;
  content += `<tr><td class="key">short name:</td><td class="value">${valueToHtml(mapSmilesToShortname(value.category))}</td></tr>`;
  content += `<tr><td class="key">smiles:</td><td class="value">${valueToHtml(value.category)}</td></tr>`;
  content += `<tr><td class="key">subset:</td><td class="value">${valueToHtml(value.subset)}</td></tr>`;
  content += `</table>`;
  content += `<div id="smiles_${value.category}" style="width:100%; height:100px; background-size: contain; background-position: center; background-repeat: no-repeat;"></div>`;

  const smiles = value.category;
  ReactionCIMEBackendFromEnv.getStructureFromSmiles(smiles, false, null).then((x) => {
    const n = document.getElementById(`smiles_${value.category}`);
    if (n != null) {
      if (x && x.length > 100) {
        // check if it is actually long enogh to be an img
        n.style.backgroundImage = `url('data:image/jpg;base64,${x}')`;
      } else {
        n.innerHTML = x;
      }
      n.title = smiles;
    }
  });
  return content;
}

function createData(feature, category, score, char) {
  return { feature, category, score, char };
}

function mapHistData(data, feature) {
  const mapped = data.map((d) => {
    return {
      feature: +d[feature],
    };
  });
  return { values: mapped };
}

function mapDensityData(allData, selectedData, feature) {
  const mappedData = allData.map((d, i) => {
    return {
      feature: +d[feature],
      selection: 'all',
    };
  });
  const mappedSelection = selectedData.map((d, i) => {
    return {
      feature: +d[feature],
      selection: 'selection',
    };
  });
  const mapped = [...mappedSelection, ...mappedData];
  return { values: mapped };
}

function mapBarChartData(allData, selectedData, feature) {
  const selectedCounts = {};
  for (let i = 0; i < selectedData.length; i++) {
    if (selectedData[i][feature] in selectedCounts) {
      selectedCounts[selectedData[i][feature]] += 1;
    } else {
      selectedCounts[selectedData[i][feature]] = 1;
    }
  }

  const allCounts = {};
  for (let i = 0; i < allData.length; i++) {
    if (allData[i][feature] in allCounts) {
      allCounts[allData[i][feature]] += 1;
    } else {
      allCounts[allData[i][feature]] = 1;
    }
  }

  const sortCountDesc = (a, b) => {
    return b.count - a.count;
  };

  const selectedBarChartData = [];
  for (const key in selectedCounts) {
    if (Object.prototype.hasOwnProperty.call(selectedCounts, key)) {
      let count = selectedCounts[key] / selectedData.length;
      count = Number.isFinite(count) ? count : 0;
      selectedBarChartData.push({ selection: 'selection', category: key, count });
    }
  }
  selectedBarChartData.sort(sortCountDesc);

  // create a map for featureCategory: id
  const categoryMap = {};
  selectedBarChartData.forEach((x, i) => {
    categoryMap[x.category] = i;
    x.id = i;
  });
  const l = selectedBarChartData.length;
  let idxCounter = l;

  const allBarChartData = [];
  for (const key in allCounts) {
    if (Object.prototype.hasOwnProperty.call(allCounts, key)) {
      let count = allCounts[key] / allData.length;
      count = Number.isFinite(count) ? count : 0;
      // apply that mapping to allBarChartData without actually having to sort it
      // make sure to check whether category in allBarChartData even exists in map, otherwise create new entry in map for new id
      let i = categoryMap[key];
      if (i == null) {
        i = idxCounter;
        categoryMap[key] = i;
        idxCounter++;
      }
      allBarChartData.push({ selection: 'all', category: key, count, id: i });
    }
  }

  const barChartData = [...allBarChartData, ...selectedBarChartData];

  return { values: barChartData };
}

const getSTD = (data) => {
  const total = data.reduce(function (a, b) {
    return a + b;
  });
  let mean = total / data.length;
  mean = Number.isFinite(mean) ? mean : 0;
  function varNumerator(value) {
    return (value - mean) * (value - mean);
  }
  let variance = data.map(varNumerator);
  variance = variance.reduce(function (a, b) {
    return a + b;
  });
  variance /= data.length;
  variance = Number.isFinite(variance) ? variance : 1;
  const std = Math.sqrt(variance);
  return std;
};

function dictionary(list) {
  const map = {};
  for (let i = 0; i < list.length; ++i) {
    for (const key in list[i]) {
      if (key in map) {
        map[key].push(list[i][key]);
      } else {
        map[key] = [list[i][key]];
      }
    }
  }
  return map;
}

function getMaxMean(data) {
  let max = Number.NEGATIVE_INFINITY;
  data = data.values;
  for (let i = 0; i < data.length; i++) {
    if (data[i].count > max) {
      max = data[i].count;
    }
  }
  return max;
}

function sortByScore(a, b) {
  if (a.score === b.score) {
    return 0;
  }

  return a.score < b.score ? 1 : -1;
}

function getProjectionColumns(legendAttributes) {
  if (legendAttributes === null) {
    return [];
  }
  const pcol = [];
  for (let i = 0; i <= legendAttributes.length; i++) {
    if (legendAttributes[i] !== undefined && legendAttributes[i].show) {
      pcol.push(legendAttributes[i].feature);
    }
  }
  return pcol;
}

function getNormalizedSTD(data, min, max) {
  if (min === max) {
    return 0;
  }
  data.forEach((x, i, self) => {
    self[i] = (+x - +min) / (+max - +min);
  });

  return getSTD(data);
}

const N_HOVER_ROWS = 5;

function genRows(vectors, aggregation, legendAttributes, dataset) {
  if (dataset === undefined) {
    return [];
  }

  const rows = [];
  const dictOfArrays = dictionary(vectors);
  const preselect = getProjectionColumns(legendAttributes);

  // loop through dict
  for (const key in dictOfArrays) {
    // filter for preselect features
    if (preselect.indexOf(key) > -1) {
      if (dataset.columns[key]?.featureType === FeatureType.Quantitative) {
        // quantitative feature
        const densityData = mapDensityData(dataset.vectors, vectors, key);
        // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
        rows.push([
          key,
          '',
          1 - getNormalizedSTD(dictOfArrays[key], dataset.columns[key].range.min, dataset.columns[key].range.max),
          <VegaDensity key={key} logLevel={vegaImport.Error} data={densityData} actions={false} tooltip={new Handler().call} />,
        ]);
      } else if (dataset.columns[key]?.featureType === FeatureType.Categorical) {
        // categorical feature
        const barData = mapBarChartData(dataset.vectors, vectors, key);
        let barChart;
        if (Object.keys(barData.values).length !== 1) {
          // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
          let tooltipOptions = {};
          if (dataset.columns[key]?.metaInformation.imgSmiles) {
            tooltipOptions = { formatTooltip: formatSMILESTooltip };
          }
          barChart = <BarChart logLevel={vegaImport.Error} data={barData} actions={false} tooltip={new Handler(tooltipOptions).call} />;
        } else {
          barChart = null;
        }
        barData.values.sort((a, b) => {
          if (a.selection === 'all') {
            return 1;
          }
          if (b.selection === 'all') {
            return -1;
          }
          return b.count - a.count;
        });
        rows.push([key, barData.values[0].category, getMaxMean(barData), barChart]);
      } else if (dataset.columns[key]?.featureType === FeatureType.Date) {
        // date feature
        const histData = mapHistData(vectors, key);
        rows.push([
          key,
          '',
          1 - getNormalizedSTD(dictOfArrays[key], dataset.columns[key].range.min, dataset.columns[key].range.max),
          <VegaDate key={key} data={histData} actions={false} tooltip={new Handler().call} />,
        ]);
      }
    }
  }

  // turn into array of dicts
  const ret = [];
  for (let i = 0; i < rows.length; i++) {
    ret.push(createData(rows[i][0], rows[i][1], rows[i][2], rows[i][3]));
  }

  // sort rows by score
  ret.sort(sortByScore);

  if (!aggregation) {
    // if it shows the hover state, we don't need to generate all rows because we can't scroll anyway -> only show top 5 charts
    return ret.slice(0, N_HOVER_ROWS);
  }

  return ret;
}

function getTable(vectors, aggregation, legendAttributes, dataset, itemLabelPlural: string) {
  const classes = makeStyles({
    table: {
      maxWidth: 288,
    },
    tableRow: {
      height: '66px',
    },
  })();
  const rows = genRows(vectors, aggregation, legendAttributes, dataset);

  return (
    <div style={{ width: '100%', maxHeight: '100%', overflowY: 'scroll' }}>
      <div
        style={{
          width: '100%',
          // overflow: "auto"
        }}
      >
        {aggregation ? (
          <Typography paddingX={2} paddingBottom={1} color="textSecondary" variant="body2" maxWidth={250}>
            The plots show distributions of feature values and are sorted by their purity.{' '}
            <Tooltip
              title={
                <Typography variant="subtitle2">
                  The plots show the value distributions of a feature overall (black outline) and the distribution of the selected {itemLabelPlural}. The
                  visualizations are sorted by the homogeneity of feature values in the subset of selected {itemLabelPlural} (i.e., measure of purity).
                </Typography>
              }
            >
              <InfoOutlined fontSize="inherit" />
            </Tooltip>
          </Typography>
        ) : (
          <Typography paddingX={1} paddingBottom={1} color="textSecondary" variant="body2" maxWidth={250}>
            Feature values of the hovered point (blue) compared to the overall value distributions (black).
          </Typography>
        )}
        <Table className={classes.table} aria-label="simple table" size="small">
          <TableHead />
          <TableBody>
            {rows.map((row) => (
              <TableRow className={classes.tableRow} key={row.feature}>
                <TableCell component="th" scope="row">
                  <div style={{ maxWidth: 200, textOverflow: 'ellipsis', overflow: 'hidden' }}>
                    {row.feature}
                    <br />
                    <b>{mapSmilesToShortname(row.category)}</b>
                  </div>
                </TableCell>
                <TableCell>{row.char}</TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
        {!aggregation && (
          <Typography paddingX={1} paddingY={1} color="textSecondary" variant="body2" maxWidth={250}>
            To show different features, adjust selection settings.
          </Typography>
        )}
      </div>
    </div>
  );
}

const mapState = (state: RootState) => {
  return {
    legendAttributes: state.genericFingerprintAttributes,
    dataset: state.dataset,
    globalLabels: state.globalLabels,
  };
};

const mapDispatch = (dispatch) => ({});

const connector = connect(mapState, mapDispatch);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  aggregate: boolean;
  selection: IVector[];
};

export const FeatureLegend = connector(({ selection, aggregate, legendAttributes, dataset, globalLabels }: Props) => {
  if (selection.length <= 0) {
    return <DefaultLegend />;
  }
  return getTable(selection, aggregate, legendAttributes, dataset, globalLabels.itemLabelPlural);
});
