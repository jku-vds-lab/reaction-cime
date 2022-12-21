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

  for (const key in showFeatures) {
    const { feature } = showFeatures[key];
    if (dataset.columns[feature]?.featureType === FeatureType.Quantitative) {
      const dataAll = await ReactionCIMEBackendFromEnv.loadDensity(dataset.info.path, feature);
      const dataLineAll = dataAll.x_vals.map((x, idx) => {
        return { feature: x, density: dataAll.y_vals[idx], selection: 'all' };
      });

      const data_hex = await ReactionCIMEBackendFromEnv.loadDensityOfHex(
        dataset.info.path,
        feature,
        workspace.xChannel,
        workspace.yChannel,
        aggregateSelection.x,
        aggregateSelection.y,
        aggregateSelection.circ_radius,
      );
      const data_line_hex = data_hex.x_vals.map((x, idx) => {
        return { feature: x, density: data_hex.y_vals[idx], selection: 'selection' };
      });

      const data_line = dataLineAll.concat(data_line_hex);
      // create vega lite chart
      var lineChart;
      if (Object.keys(data_line).length > 1) {
        // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
        lineChart = <AreaChart logLevel={vegaImport.Error} data={{ values: data_line }} actions={false} />;
      } else {
        lineChart = null;
      }
      rows.push({ feature, category: '', char: lineChart, score: 1 - data_hex.norm_sd });
    } else {
      const { bar_data, cat_list, score } = await handleCategoricalData(dataset, feature, workspace, aggregateSelection);

      // create vega lite chart
      var barChart;
      if (Object.keys(bar_data).length > 1) {
        // logLevel={vegaImport.Debug} | {vegaImport.Warn} | {vegaImport.Error} | {vegaImport.None} | {vegaImport.Info}
        let tooltip_options = {};
        if (dataset.columns[feature]?.metaInformation.imgSmiles) {
          tooltip_options = { formatTooltip: formatSMILESTooltip };
        }
        barChart = <BarChart logLevel={vegaImport.Error} data={{ values: bar_data }} actions={false} tooltip={new Handler(tooltip_options).call} />;
      } else {
        barChart = null;
      }
      rows.push({ feature, category: cat_list[0], char: barChart, score });
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
  const data_all = await ReactionCIMEBackendFromEnv.loadCategoryCount(dataset.info.path, feature);
  const sum_all = data_all.reduce(function (acc, row) {
    return acc + row.count;
  }, 0);
  const bar_data_all = data_all.map((row) => {
    return { selection: 'all', category: row[feature], count: row.count / sum_all };
  });

  const data_hex = await ReactionCIMEBackendFromEnv.loadCategoryCountOfHex(
    dataset.info.path,
    feature,
    workspace.xChannel,
    workspace.yChannel,
    aggregateSelection.x,
    aggregateSelection.y,
    aggregateSelection.circ_radius,
  );
  const sum_hex = data_hex.reduce(function (acc, row) {
    return acc + row.count;
  }, 0);
  const bar_data_hex = data_hex.map((row) => {
    return { selection: 'selection', category: row[feature], count: row.count / sum_hex };
  });

  // sort cat_list by count
  const cat_list = bar_data_all.map((row) => row.category);
  cat_list.sort((a, b) => {
    const bar_hex_a = bar_data_hex.find((row) => row.category === a);
    const bar_hex_b = bar_data_hex.find((row) => row.category === b);
    if (bar_hex_a && bar_hex_b) {
      // first sort by aggregate count
      return bar_hex_b.count - bar_hex_a.count;
    }
    if (bar_hex_a) {
      // if there exists this value
      return -1;
    }
    if (bar_hex_b) {
      // if there exists this value
      return 1;
    }

    // if feature not present in aggregated data
    const bar_all_a = bar_data_all.find((row) => row.category === a);
    const bar_all_b = bar_data_all.find((row) => row.category === b);
    return bar_all_b.count - bar_all_a.count;
  });

  // add ids to the data; i.e. instances with the same feature values must have same id
  const bar_data = [];
  for (const i in cat_list) {
    const cat = cat_list[i];

    const bar_all = bar_data_all.find((row) => row.category === cat);
    if (bar_all != null) {
      bar_data.push({ ...bar_all, id: i });
    }
    const bar_hex = bar_data_hex.find((row) => row.category === cat);
    if (bar_hex != null) {
      bar_data.push({ ...bar_hex, id: i });
    }
  }

  const score = bar_data_hex.find((row) => row.category === cat_list[0]).count; // maximum bar count is score
  return { bar_data, cat_list, score };
}
