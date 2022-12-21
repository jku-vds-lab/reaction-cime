import { Table, TableBody, TableCell, TableHead, TableRow } from '@mui/material';
import { makeStyles } from '@mui/styles';
import { Dataset, FeatureType, IProjection } from 'projection-space-explorer';
import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import * as vegaImport from 'vega';
import { Handler } from 'vega-tooltip';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { AppState } from '../../State/Store';
import BarChart from './VegaHelpers/BarChart';
import AreaChart from './VegaHelpers/AreaChart';
import { formatSMILESTooltip } from './FeatureLegend';
import { mapSmilesToShortname } from '../../Utility/Utils';

async function genRows(
  aggregateSelection: { x: number; y: number; circ_radius: number },
  aggregation: boolean,
  legendAttributes: { feature: string; show: boolean }[],
  dataset: Dataset,
  workspace: IProjection,
  setRows: any,
) {
  const showFeatures = legendAttributes.filter((row) => row.show);
  const rows = [];

  // eslint-disable-next-line guard-for-in
  for (const key in showFeatures) {
    const { feature } = showFeatures[key];
    if (dataset.columns[feature]?.featureType === FeatureType.Quantitative) {
      const dataAll = await ReactionCIMEBackendFromEnv.loadDensity(dataset.info.path, feature);
      const dataLineAll = dataAll.x_vals.map((x, idx) => {
        return { feature: x, density: dataAll.y_vals[idx], selection: 'all' };
      });

      const dataHex = await ReactionCIMEBackendFromEnv.loadDensityOfHex(
        dataset.info.path,
        feature,
        workspace.xChannel,
        workspace.yChannel,
        aggregateSelection.x,
        aggregateSelection.y,
        aggregateSelection.circ_radius,
      );
      const dataLineHex = dataHex.x_vals.map((x, idx) => {
        return { feature: x, density: dataHex.y_vals[idx], selection: 'selection' };
      });

      const dataLine = dataLineAll.concat(dataLineHex);
      // create vega lite chart
      let lineChart;
      if (Object.keys(dataLine).length > 1) {
        // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
        lineChart = <AreaChart logLevel={vegaImport.Error} data={{ values: dataLine }} actions={false} />;
      } else {
        lineChart = null;
      }
      rows.push({ feature, category: '', char: lineChart, score: 1 - dataHex.norm_sd });
    } else {
      const { bar_data: barData, cat_list: catList, score } = await handleCategoricalData(dataset, feature, workspace, aggregateSelection);

      // create vega lite chart
      let barChart;
      if (Object.keys(barData).length > 1) {
        // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
        let tooltipOptions = {};
        if (dataset.columns[feature]?.metaInformation.imgSmiles) {
          tooltipOptions = { formatTooltip: formatSMILESTooltip };
        }
        barChart = <BarChart logLevel={vegaImport.Error} data={{ values: barData }} actions={false} tooltip={new Handler(tooltipOptions).call} />;
      } else {
        barChart = null;
      }
      rows.push({ feature, category: catList[0], char: barChart, score });
    }
  }

  // sort rows by score
  rows.sort((a, b) => b.score - a.score);

  setRows(rows);
}

const mapStateToProps = (state: AppState) => ({
  // aggregateSelection: state.selection?.currentAggregateSelection,
  legendAttributes: state.genericFingerprintAttributes,
  dataset: state.dataset,
  workspace: state.multiples.multiples.entities[state.multiples.active].attributes.workspace,
});
const mapDispatchToProps = (dispatch: any) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  aggregate: boolean;
  aggregateSelection: any;
};

const useStyles = makeStyles({
  table: {
    maxWidth: 288,
  },
  tableRow: {
    height: '66px',
  },
});

export const AggregateLegend = connector(({ aggregate, aggregateSelection, legendAttributes, dataset, workspace }: Props) => {
  const classes = useStyles();
  const [rows, setRows] = React.useState([]);

  React.useEffect(() => {
    genRows(aggregateSelection, aggregate, legendAttributes, dataset, workspace as IProjection, setRows);
    // only need to update if we have a new aggregate selection
    // eslint-disable-next-line
  }, [aggregateSelection]);
  return (
    <div style={{ width: '100%', maxHeight: '100%', overflowY: 'scroll' }}>
      <div
        style={{
          width: '100%',
          // overflow: "auto"
        }}
      >
        <Table className={classes.table} aria-label="simple table" size="small">
          <TableHead />
          <TableBody>
            {rows.map((row) => (
              <TableRow className={classes.tableRow} key={row.feature}>
                <TableCell component="th" scope="row">
                  <div style={{ maxWidth: 200 }}>
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
      </div>
    </div>
  );
});

async function handleCategoricalData(
  dataset: Dataset,
  feature: string,
  workspace: IProjection,
  aggregateSelection: { x: number; y: number; circ_radius: number },
) {
  const dataAll = await ReactionCIMEBackendFromEnv.loadCategoryCount(dataset.info.path, feature);
  const sumAll = dataAll.reduce(function (acc, row) {
    return acc + row.count;
  }, 0);
  const barDataAll = dataAll.map((row) => {
    return { selection: 'all', category: row[feature], count: row.count / sumAll };
  });

  const dataHex = await ReactionCIMEBackendFromEnv.loadCategoryCountOfHex(
    dataset.info.path,
    feature,
    workspace.xChannel,
    workspace.yChannel,
    aggregateSelection.x,
    aggregateSelection.y,
    aggregateSelection.circ_radius,
  );
  const sumHex = dataHex.reduce(function (acc, row) {
    return acc + row.count;
  }, 0);
  const barDataHex = dataHex.map((row) => {
    return { selection: 'selection', category: row[feature], count: row.count / sumHex };
  });

  // sort cat_list by count
  const catList = barDataAll.map((row) => row.category);
  catList.sort((a, b) => {
    const barHexA = barDataHex.find((row) => row.category === a);
    const barHexB = barDataHex.find((row) => row.category === b);
    if (barHexA && barHexB) {
      // first sort by aggregate count
      return barHexB.count - barHexA.count;
    }
    if (barHexA) {
      // if there exists this value
      return -1;
    }
    if (barHexB) {
      // if there exists this value
      return 1;
    }

    // if feature not present in aggregated data
    const barAllA = barDataAll.find((row) => row.category === a);
    const barAllB = barDataAll.find((row) => row.category === b);
    return barAllB.count - barAllA.count;
  });

  // add ids to the data; i.e. instances with the same feature values must have same id
  const barData = [];
  catList.forEach((cat, i) => {
    const barAll = barDataAll.find((row) => row.category === cat);
    if (barAll != null) {
      barData.push({ ...barAll, id: i });
    }
    const barHex = barDataHex.find((row) => row.category === cat);
    if (barHex != null) {
      barData.push({ ...barHex, id: i });
    }
  });

  const score = barDataHex.find((row) => row.category === catList[0]).count; // maximum bar count is score
  return { bar_data: barData, cat_list: catList, score };
}
