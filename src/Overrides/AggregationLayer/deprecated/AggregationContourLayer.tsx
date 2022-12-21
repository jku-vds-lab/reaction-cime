import { useCancellablePromise } from 'projection-space-explorer';
import * as React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import * as THREE from 'three';
import * as _ from 'lodash';
import * as vsup from 'vsup';
import * as d3 from 'd3v5';
import { ReactionCIMEBackendFromEnv } from '../../../Backend/ReactionCIMEBackend';
import { AppState } from '../../../State/Store';
import { AggregateDataset } from '../AggregateDataset';
import { LoadingIndicatorDialog } from '../../Dataset/LoadingIndicatorDialog';
import { GLContour } from './GLContour';
import { AggregateActions } from '../../../State/AggregateSettingsDuck';

const retrieveColorscale = (
  aggregateDataset: AggregateDataset,
  value_col: string,
  uncertainty_col: string,
  setAggregateColorMapScale: any,
  aggregateSettings: any,
) => {
  if (!Object.keys(aggregateDataset.columns).includes(value_col)) {
    return [null, null];
  }

  // visit https://github.com/uwdata/vsup for more info about VSUP vs bivariate colorscale encoding
  const vDom = [aggregateDataset.columns[value_col].range.min, aggregateDataset.columns[value_col].range.max] as [number, number];

  let scale;
  let quantization;
  if (aggregateDataset.columns[uncertainty_col] != null) {
    const uDom = [aggregateDataset.columns[uncertainty_col].range.min, aggregateDataset.columns[uncertainty_col].range.max];

    if (aggregateSettings?.useVSUP) {
      quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
    } else {
      quantization = vsup.squareQuantization(4).valueDomain(vDom).uncertaintyDomain(uDom);
    }
    scale = vsup.scale().quantize(quantization).range(d3[aggregateSettings?.colorscale]);
  } else {
    scale = d3.scaleQuantize().domain(vDom).range(d3.quantize(d3[aggregateSettings.colorscale], 8));
  }
  setAggregateColorMapScale(scale);
  return scale;
};

const createContours = (dataset, value_col, scale) => {
  if (!Object.keys(dataset.columns).includes(value_col)) {
    return [null, null];
  }
  // if(aggregateDataset.columns[uncertainty_col] != null){
  //     var uncertainty_arr = dataset.vectors.map((row) => row[uncertainty_col]);
  // }

  const valueArr = dataset.vectors.map((row) => row[value_col]);
  // console.log(value_arr)
  // let [width, height, x, y] = [100,100,0,0];
  // if(dataset.bounds){
  //     width = dataset.bounds.x["max"]-dataset.bounds.x["min"];
  //     height = dataset.bounds.y["max"]-dataset.bounds.y["min"];

  //     // need to set x and y because the center is not at 0 0 if max and min are unequal
  //     x = (dataset.bounds.x["max"]+dataset.bounds.x["min"])/2;
  //     y = (dataset.bounds.y["max"]+dataset.bounds.y["min"])/2;
  // }

  const xAxis = d3
    .scaleLinear()
    .range([0, Math.sqrt(valueArr.length)])
    .domain([dataset.bounds.x.min, dataset.bounds.x.max]);

  const yAxis = d3
    .scaleLinear()
    .range([0, Math.sqrt(valueArr.length)])
    .domain([dataset.bounds.y.min, dataset.bounds.y.max]);

  const min = Math.min(...valueArr);
  const max = Math.max(...valueArr);
  console.log(min, max);
  const noNanArr = valueArr.map((val) => (Number.isNaN(parseFloat(val)) ? -Number.MAX_VALUE : val));
  console.log(noNanArr);
  const contours = d3
    .contours()
    // .thresholds(11)
    .thresholds(d3.range(min, max, (max - min) / 11))
    .smooth(true)
    .size([Math.sqrt(valueArr.length), Math.sqrt(valueArr.length)])(noNanArr);

  const lines = [];
  contours.forEach((contour) => {
    // const coordinates = contour.coordinates[0][0]
    const coordinatesList = contour.coordinates;
    const { value } = contour;
    // let material = new THREE.LineBasicMaterial({ color: scale(value) })
    const material = new THREE.MeshBasicMaterial({ color: scale(value), side: THREE.DoubleSide });
    coordinatesList.forEach((arr) => {
      const coordinates = arr[0];
      const shape = new THREE.Shape();
      shape.moveTo(xAxis.invert(coordinates[0][0]), yAxis.invert(coordinates[0][1]));
      for (let i = 0; i < coordinates.length; i++) {
        const cur = coordinates[i];
        shape.lineTo(xAxis.invert(cur[0]), yAxis.invert(cur[1]));
      }

      const geo = new THREE.ShapeGeometry(shape);
      const line = new THREE.Mesh(geo, material);

      line.visible = true;
      lines.push(line);
    });
  });

  return lines;
};

const mapStateToProps = (state: AppState) => ({
  // aggregateColor: state.aggregateSettings?.aggregateColor,
  aggregateColor: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.aggregateSettings?.colormapSettings.aggregateColor,
  poiDataset: state.dataset,
  viewTransform: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.viewTransform,
  // aggregateSettings: state.aggregateSettings,
  aggregateSettings: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.aggregateSettings,
});
const mapDispatchToProps = (dispatch: any) => ({
  // setAggregateDataset: dataset => dispatch(setAggregateDatasetAction(dataset)),
  setAggregateColorMapScale: (legend) => dispatch(AggregateActions.setAggregateColorMapScale(legend)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux;

const loadingArea = 'global_loading_indicator_aggregation_ds';
export const AggregationContourLayer = connector(
  ({ aggregateColor, poiDataset, viewTransform, setAggregateColorMapScale, aggregateSettings }: AggregationLayerProps) => {
    if (poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === 'None') {
      return null;
    }

    const [lines, setLines] = React.useState(null);

    const [aggregateDataset, setAggregateDataset] = React.useState(null);

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    React.useEffect(() => {
      // setAggregateDataset(null)
      const abortController = new AbortController(); // TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
      // load the basic aggregateDataset with the high-level overview information
      ReactionCIMEBackendFromEnv.loadAggCSV(
        (dataset) => {
          setAggregateDataset(new AggregateDataset(dataset));
        },
        poiDataset.info.path,
        aggregateColor.value_col,
        aggregateColor.uncertainty_col,
        aggregateColor.cache_cols,
        0,
        null,
        cancellablePromise,
        abortController,
        loadingArea,
      );
    }, [aggregateColor, poiDataset.info.path]);

    React.useEffect(() => {
      if (aggregateDataset && aggregateDataset.vectors) {
        retrieveColorscale(aggregateDataset, aggregateColor.value_col, aggregateColor.uncertainty_col, setAggregateColorMapScale, aggregateSettings);
      }
    }, [aggregateDataset, aggregateColor, aggregateSettings?.colormapSettings.colorscale, aggregateSettings?.colormapSettings.useVSUP]);

    React.useEffect(() => {
      if (aggregateSettings?.colormapSettings.scale_obj != null && aggregateDataset && aggregateDataset.vectors) {
        const lines = createContours(aggregateDataset, aggregateColor.value_col, aggregateSettings?.colormapSettings.scale_obj);
        setLines(lines);
      } else {
        setLines(null);
      }
    }, [aggregateSettings?.colormapSettings.scale_obj, aggregateSettings?.colormapSettings.valueFilter]);

    return (
      <div>
        <LoadingIndicatorDialog
          handleClose={() => {
            cancelPromises();
          }}
          area={loadingArea}
        />
        {lines && lines.length > 0 && <GLContour lines={lines} />}
      </div>
    );
  },
);
