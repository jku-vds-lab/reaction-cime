import { useCancellablePromise } from 'projection-space-explorer';
import * as React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import * as THREE from 'three';
import * as _ from 'lodash';
import { ReactionCIMEBackendFromEnv } from '../../../Backend/ReactionCIMEBackend';
import { AppState } from '../../../State/Store';
import { AggregateDataset } from '../AggregateDataset';
import { LoadingIndicatorDialog } from '../../Dataset/LoadingIndicatorDialog';
import { GLHeatmap } from './GLHeatmap';
import { AggregateActions } from '../../../State/AggregateSettingsDuck';
import { convertToRgb } from '../../../Utility/Utils';

const retrieveInformationFromAggDataset = (
  aggregateDataset: AggregateDataset,
  dataset: AggregateDataset,
  value_col: string,
  uncertainty_col: string,
  valueFilter: string[],
  scale: any,
) => {
  let uncertaintyArr: any[];
  if (aggregateDataset.columns[uncertainty_col] != null) {
    uncertaintyArr = dataset.vectors.map((row) => row[uncertainty_col]);
  }

  const valueArr = dataset.vectors.map((row) => row[value_col]);
  const bgRGBA = new Uint8Array(valueArr.length * 4);
  for (let i = 0; i < valueArr.length; i++) {
    // set opacity to 0 if value is not given
    if (valueArr[i] === undefined || Number.isNaN(valueArr[i]) || valueArr[i] === '') {
      bgRGBA[4 * i + 3] = 0;
    } else {
      let color;
      if (uncertaintyArr) {
        color = scale(valueArr[i], uncertaintyArr[i]);
      } else {
        color = scale(valueArr[i], 0);
      }

      if (valueFilter != null && valueFilter.length > 0 && !valueFilter.includes(color)) {
        // if there is a filter, and the current value is not within, we set it to transparent
        bgRGBA[4 * i + 3] = 0;
      } else {
        const rgb = convertToRgb(color);

        // RGB from 0 to 255
        bgRGBA[4 * i] = rgb.r;
        bgRGBA[4 * i + 1] = rgb.g;
        bgRGBA[4 * i + 2] = rgb.b;

        // OPACITY
        bgRGBA[4 * i + 3] = 255;
      }
    }
  }
  const bgDataTex = new THREE.DataTexture(bgRGBA, Math.sqrt(valueArr.length), Math.sqrt(valueArr.length), THREE.RGBAFormat);
  // bgDataTex.magFilter = THREE.LinearMipMapLinearFilter; // this causes border artefacts
  bgDataTex.magFilter = THREE.NearestFilter; // this makes it discrete
  // bgDataTex.magFilter = THREE.LinearFilter; // this causes border artefacts
  // bgDataTex.minFilter = THREE.LinearMipMapLinearFilter; // this causes everything to go black
  // bgDataTex.minFilter = THREE.LinearFilter;
  bgDataTex.minFilter = THREE.NearestFilter;

  let [width, height, x, y] = [100, 100, 0, 0];
  if (dataset.bounds) {
    width = dataset.bounds.x.max - dataset.bounds.x.min;
    height = dataset.bounds.y.max - dataset.bounds.y.min;

    // need to set x and y because the center is not at 0 0 if max and min are unequal
    x = (dataset.bounds.x.max + dataset.bounds.x.min) / 2;
    y = (dataset.bounds.y.max + dataset.bounds.y.min) / 2;
  }

  return [bgDataTex, { width, height, x, y }];
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
  setValueRange: (range) => dispatch(AggregateActions.setValueRange(range)),
  setUncertaintyRange: (range) => dispatch(AggregateActions.setUncertaintyRange(range)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux;

const loadingArea = 'global_loading_indicator_aggregation_ds';
const AggregationLayer = connector(
  ({ aggregateColor, poiDataset, viewTransform, aggregateSettings, setValueRange, setUncertaintyRange }: AggregationLayerProps) => {
    if (poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === 'None') {
      return null;
    }

    const [textures, setTextures] = React.useState(null);
    const [sizes, setSizes] = React.useState(null);

    const [aggregateDataset, setAggregateDataset] = React.useState(null);
    const [aggregateDatasetZoomed, setAggregateDatasetZoomed] = React.useState(null);

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    // eslint-disable-next-line react-hooks/exhaustive-deps
    const debouncedLoadAggDataset = React.useCallback(
      _.debounce(
        (newViewTransform) => {
          cancelPromises();
          const abortController = new AbortController(); // TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?

          const range = {
            x: {
              min: Math.round(-newViewTransform.width / newViewTransform.zoom / 2 + newViewTransform.centerX),
              max: Math.round(newViewTransform.width / newViewTransform.zoom / 2 + newViewTransform.centerX),
            },
            y: {
              min: Math.round(-newViewTransform.height / newViewTransform.zoom / 2 + newViewTransform.centerY),
              max: Math.round(newViewTransform.height / newViewTransform.zoom / 2 + newViewTransform.centerY),
            },
          };

          // take 150% of the boundaries, such that we have clean borders
          range.x.min -= Math.abs(range.x.min / 2);
          range.x.max += Math.abs(range.x.max / 2);
          range.y.min -= Math.abs(range.y.min / 2);
          range.y.max += Math.abs(range.y.max / 2);

          // load zoomed version of the aggregate dataset
          ReactionCIMEBackendFromEnv.loadAggCSV(
            (dataset) => {
              setAggregateDatasetZoomed(new AggregateDataset(dataset));
            },
            poiDataset.info.path,
            aggregateColor.value_col,
            aggregateColor.uncertainty_col,
            aggregateColor.cache_cols,
            0,
            range,
            cancellablePromise,
            abortController,
            'None',
          );
        },
        500,
        { leading: false, trailing: true },
      ), // leading: execute function at the beginning of the events; trailing: execute function at the end of the events; maxWait: maximum time the function is allowed to be delayed
      [aggregateColor],
    );

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

      // reset the zoomed version of the dataset
      setAggregateDatasetZoomed(null);
      // eslint-disable-next-line
    }, [aggregateColor, poiDataset.info.path]);

    React.useEffect(() => {
      if (aggregateDataset != null)
        // only start loading, if we already have a base dataset
        debouncedLoadAggDataset(viewTransform);
      // eslint-disable-next-line
    }, [aggregateDataset, viewTransform]);

    React.useEffect(() => {
      if (aggregateDataset && aggregateDataset.vectors && aggregateSettings?.advancedSettings.deriveRange) {
        if (Object.keys(aggregateDataset.columns).includes(aggregateColor.value_col)) {
          setValueRange(aggregateDataset.columns[aggregateColor.value_col].range);

          if (aggregateDataset.columns[aggregateColor.uncertainty_col] != null) {
            setUncertaintyRange(aggregateDataset.columns[aggregateColor.uncertainty_col].range);
          } else {
            setUncertaintyRange(null);
          }
        }
      }
      // aggregateColor --> aggregateColor has direct influence on aggregateDataset through "setAggregateDataset" in the useEffect above
      // eslint-disable-next-line
    }, [aggregateDataset, aggregateSettings?.advancedSettings.deriveRange]);

    React.useEffect(() => {
      if (aggregateSettings.colormapSettings.scale_obj != null) {
        if (aggregateDataset && aggregateDataset.vectors) {
          if (Object.keys(aggregateDataset.columns).includes(aggregateColor.value_col)) {
            const newSizes = [null, null];
            const newTextures = [null, null];
            const [texture, size] = retrieveInformationFromAggDataset(
              aggregateDataset,
              aggregateDataset,
              aggregateColor.value_col,
              aggregateColor.uncertainty_col,
              aggregateSettings?.colormapSettings.valueFilter,
              aggregateSettings.colormapSettings.scale_obj,
            );
            newSizes[0] = size;
            newTextures[0] = texture;

            if (aggregateDatasetZoomed && aggregateDatasetZoomed.vectors) {
              const retrieved = retrieveInformationFromAggDataset(
                aggregateDataset,
                aggregateDatasetZoomed,
                aggregateColor.value_col,
                aggregateColor.uncertainty_col,
                aggregateSettings?.colormapSettings.valueFilter,
                aggregateSettings.colormapSettings.scale_obj,
              );
              newSizes[1] = retrieved[1];
              newTextures[1] = retrieved[0];
            }
            setSizes(newSizes);
            setTextures(newTextures);
          }
        } else {
          setSizes(null);
          setTextures(null);
        }
      }

      // scale_obj is directly dependent on uncertainty_col and value_col
      // eslint-disable-next-line
    }, [aggregateSettings?.colormapSettings.scale_obj, aggregateDataset, aggregateDatasetZoomed, aggregateSettings?.colormapSettings.valueFilter]);

    return (
      <div>
        <LoadingIndicatorDialog
          handleClose={() => {
            cancelPromises();
          }}
          area={loadingArea}
        />
        {textures && textures.length > 0 && sizes && sizes.length > 0 && (
          <GLHeatmap
            // texture={new THREE.TextureLoader().load('https://threejsfundamentals.org/threejs/resources/images/wall.jpg')}
            textures={textures}
            sizes={sizes}
          />
        )}
      </div>
    );
  },
);

export { AggregationLayer };
