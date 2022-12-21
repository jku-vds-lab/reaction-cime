import { IProjection, useCancellablePromise } from 'projection-space-explorer';
import * as React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import * as THREE from 'three';
import * as _ from 'lodash';
import { EntityId } from '@reduxjs/toolkit';
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { AppState } from '../../State/Store';
import { AggregateDataset } from './AggregateDataset';
import { LoadingIndicatorDialog } from '../Dataset/LoadingIndicatorDialog';
import { AggregateActions } from '../../State/AggregateSettingsDuck';
import { GLHexagons } from './GLHexagons';
import { PSE_BLUE } from '../../Utility/Utils';
import { setCurrentAggregateSelection } from '../../State/SelectionDuck';
import { ReactionVector } from '../../State/interfaces';

const createHexagons = (
  dataset: AggregateDataset,
  xChannel: string,
  yChannel: string,
  value_col: string,
  uncertainty_col: string,
  scale_obj: any,
  valueFilter: string[],
) => {
  const hexagons = [];

  dataset.vectors.forEach((row: ReactionVector) => {
    if (row.hex === 'False') {
      console.log('wrong point');
      const materialWrong = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.DoubleSide });
      const geometryWrong = new THREE.CircleGeometry(1, 100);
      const objectWrong = new THREE.Mesh(geometryWrong, materialWrong);
      objectWrong.position.x = row[xChannel];
      objectWrong.position.y = row[yChannel];

      hexagons.push(objectWrong);
    } else if (row[value_col] === undefined || Number.isNaN(row[value_col]) || row[value_col] === '') {
      console.warn('empty block');
    } else {
      let color;
      if (dataset.columns[uncertainty_col] != null) {
        color = scale_obj(row[value_col], row[uncertainty_col]);
      } else {
        color = scale_obj(row[value_col], 0);
      }
      if (valueFilter == null || valueFilter.length <= 0 || valueFilter.includes(color)) {
        // if there is a filter, and the current value is not within, we ignore this hexagon
        const material = new THREE.MeshBasicMaterial({ color, side: THREE.DoubleSide });
        const geometry = new THREE.CircleGeometry(row.circ_radius - row.circ_radius * 0.04, 6);
        const object = new THREE.Mesh(geometry, material);
        object.position.x = row[xChannel];
        object.position.y = row[yChannel];

        hexagons.push(object);
      }
    }
  });

  return hexagons;
};

// test if point is inside hexagon adapted from http://www.playchilla.com/how-to-check-if-a-point-is-inside-a-hexagon
function isInsideHex(point_x: number, point_y: number, hex_x: number, hex_y: number, radius: number, circ_radius: number) {
  const pointq2x = Math.abs(point_x - hex_x); // transform the test point locally and to quadrant 2
  const pointq2y = Math.abs(point_y - hex_y); // transform the test point locally and to quadrant 2
  if (pointq2x > circ_radius || pointq2y > radius) {
    return false;
  }
  return radius * circ_radius - radius * pointq2x - (circ_radius / 2) * pointq2y >= 0; // finally the dot product can be reduced to this due to the hexagon symmetry
}

// test if point is outside of a bounding box, including a certain threshold
function isOutsideBoundingBox(position: { x: number; y: number }, dataset: AggregateDataset, threshold: number, xChannel: string, yChannel: string) {
  if (!(xChannel in dataset.columns) || !(yChannel in dataset.columns)) return false;
  return (
    position.x < dataset.columns[xChannel].range.min - threshold ||
    position.x > dataset.columns[xChannel].range.max + threshold ||
    position.y < dataset.columns[yChannel].range.min - threshold ||
    position.y > dataset.columns[yChannel].range.max + threshold
  );
}

const mapStateToProps = (state: AppState) => ({
  // aggregateColor: state.aggregateSettings?.aggregateColor,
  // aggregateColor: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.aggregateSettings?.aggregateColor,
  poiDataset: state.dataset,
  smallMultiples: state.multiples.multiples.entities,
  // viewTransform: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.viewTransform,
  // aggregateSettings: state.aggregateSettings,
  // aggregateSettings: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.aggregateSettings,
  mouseMove: state.mouseInteractionHooks?.mousemove,
  mouseClick: state.mouseInteractionHooks?.mouseclick,
});
const mapDispatchToProps = (dispatch: any) => ({
  setValueRange: (range) => dispatch(AggregateActions.setValueRange(range)),
  setUncertaintyRange: (range) => dispatch(AggregateActions.setUncertaintyRange(range)),
  setCurrentAggregateSelectionFn: (selection) => dispatch(setCurrentAggregateSelection(selection)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
  multipleId: EntityId;
};

const loadingArea = 'global_loading_indicator_aggregation_ds';
export const HexAggregationLayer = connector(
  ({
    setCurrentAggregateSelectionFn,
    poiDataset,
    smallMultiples,
    setValueRange,
    setUncertaintyRange,
    mouseMove,
    mouseClick,
    multipleId,
  }: AggregationLayerProps) => {
    const { aggregateSettings } = smallMultiples[multipleId].attributes;
    const aggregateColor = aggregateSettings?.colormapSettings.aggregateColor;
    const { viewTransform } = smallMultiples[multipleId].attributes;
    const workspace = smallMultiples[multipleId].attributes.workspace as IProjection;
    const xChannel = workspace.xChannel == null ? 'x' : workspace.xChannel;
    const yChannel = workspace.yChannel == null ? 'y' : workspace.yChannel;

    if (poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === 'None') {
      return null;
    }

    const [hexagons, setHexagons] = React.useState(null);
    const [hoverElement, setHoverElement] = React.useState(null);
    const [selectElement, setSelectElement] = React.useState(null);
    const [aggregateDataset, setAggregateDataset] = React.useState<AggregateDataset>(null);
    const [datasetValueRange, setDatasetValueRange] = React.useState(null);
    const [datasetUncertaintyRange, setDatasetUncertaintyRange] = React.useState(null);

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    React.useEffect(() => {
      if (aggregateColor?.value_col != null) {
        ReactionCIMEBackendFromEnv.loadValueRange(poiDataset.info.path, aggregateColor.value_col).then((response) => {
          setDatasetValueRange(response);
        });
        if (aggregateColor?.uncertainty_col != null) {
          ReactionCIMEBackendFromEnv.loadValueRange(poiDataset.info.path, aggregateColor.uncertainty_col).then((response) => {
            setDatasetUncertaintyRange(response);
          });
        } else {
          setDatasetUncertaintyRange(null);
        }
      }
    }, [aggregateColor, poiDataset.info.path]);

    React.useEffect(() => {
      if (aggregateSettings?.advancedSettings.deriveRange) {
        if (datasetValueRange != null) {
          setValueRange(datasetValueRange);
          setUncertaintyRange(datasetUncertaintyRange);
        }
      }
      // eslint-disable-next-line
    }, [datasetValueRange, datasetUncertaintyRange, aggregateSettings?.advancedSettings.deriveRange]);

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

          // load the aggregateDataset
          ReactionCIMEBackendFromEnv.loadHexAgg(
            (dataset) => {
              setAggregateDataset(new AggregateDataset(dataset, xChannel, yChannel));
            },
            poiDataset.info.path,
            xChannel,
            yChannel,
            aggregateColor.value_col,
            aggregateColor.uncertainty_col,
            aggregateColor.cache_cols,
            0,
            aggregateSettings?.advancedSettings.aggregationMethod,
            range,
            cancellablePromise,
            abortController,
            loadingArea,
          );
        },
        500,
        { leading: false, trailing: true },
      ), // leading: execute function at the beginning of the events; trailing: execute function at the end of the events; maxWait: maximum time the function is allowed to be delayed
      [poiDataset?.info?.path, aggregateColor, aggregateSettings, xChannel, yChannel],
    );

    React.useEffect(() => {
      // setAggregateDataset(null)
      debouncedLoadAggDataset(viewTransform);
      // eslint-disable-next-line
    }, [aggregateColor, poiDataset.info.path, aggregateSettings?.advancedSettings.aggregationMethod, viewTransform, xChannel, yChannel]);

    React.useEffect(() => {
      if (aggregateSettings?.colormapSettings.scale_obj != null) {
        if (aggregateDataset && aggregateDataset.vectors) {
          if (Object.keys(aggregateDataset.columns).includes(aggregateColor.value_col)) {
            const hexs = createHexagons(
              aggregateDataset,
              xChannel,
              yChannel,
              aggregateColor.value_col,
              aggregateColor.uncertainty_col,
              aggregateSettings?.colormapSettings.scale_obj,
              aggregateSettings?.colormapSettings.valueFilter,
            );
            setHexagons(hexs);
          }
        } else {
          setHexagons(null);
        }
      }

      // aggregateColor --> aggregateColor has direct influence on aggregateDataset through "setAggregateDataset" in the useEffect above
      // eslint-disable-next-line
    }, [aggregateSettings?.colormapSettings.scale_obj, aggregateDataset, aggregateSettings?.colormapSettings.valueFilter]);

    React.useEffect(() => {
      if (mouseMove != null && !mouseMove.event_used && aggregateDataset != null) {
        const circRadius = aggregateDataset.columns.circ_radius.range.max;
        // check if it is outside of the aggregation bounding box
        if (isOutsideBoundingBox(mouseMove, aggregateDataset, circRadius, xChannel, yChannel)) {
          setHoverElement(null);
        } else {
          // if it is inside, we need to check each hexagon, if coordinates are inside
          const radius = (Math.sqrt(3) / 2) * circRadius; // https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
          let foundHex = false;
          aggregateDataset.vectors.forEach((row: ReactionVector) => {
            const isInside = isInsideHex(mouseMove.x, mouseMove.y, row[xChannel], row[yChannel], radius, circRadius);
            if (isInside) {
              foundHex = true;
              if (hoverElement == null || hoverElement.position.x !== row[xChannel] || hoverElement.position.y !== row[yChannel]) {
                const material = new THREE.MeshBasicMaterial({ color: PSE_BLUE, side: THREE.DoubleSide, transparent: true, opacity: 1.0 });
                const geometry = new THREE.CircleGeometry(row.circ_radius - row.circ_radius * 0.04, 6);
                const object = new THREE.Mesh(geometry, material);
                object.position.x = row[xChannel];
                object.position.y = row[yChannel];
                setHoverElement(object);
              }
            }
          });

          if (!foundHex) {
            // if position is not inside a hexagon, we set hover to null
            setHoverElement(null);
          }
        }
      } else {
        setHoverElement(null);
      }

      // eslint-disable-next-line
    }, [aggregateDataset, mouseMove]);

    React.useEffect(() => {
      if (mouseClick != null && !mouseClick.event_used && aggregateDataset != null) {
        const circRadius = aggregateDataset.columns.circ_radius.range.max;
        // check if it is outside of the aggregation bounding box
        if (isOutsideBoundingBox(mouseClick, aggregateDataset, circRadius, xChannel, yChannel)) {
          setSelectElement(null);
          setCurrentAggregateSelectionFn(null);
        } else {
          // if it is inside, we need to check each hexagon, if coordinates are inside
          const radius = (Math.sqrt(3) / 2) * circRadius; // https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
          let foundHex = false;
          aggregateDataset.vectors.forEach((row: ReactionVector) => {
            const isInside = isInsideHex(mouseClick.x, mouseClick.y, row[xChannel], row[yChannel], radius, circRadius);
            if (isInside) {
              foundHex = true;
              if (selectElement == null || selectElement.position.x !== row.x || selectElement.position.y !== row.y) {
                const material = new THREE.MeshBasicMaterial({ color: PSE_BLUE, side: THREE.DoubleSide, transparent: true, opacity: 1.0 });
                const tempRadius = Number(row.circ_radius) + Number(row.circ_radius) * 0.04;
                const geometry = new THREE.CircleGeometry(tempRadius, 6);
                const object = new THREE.Mesh(geometry, material);
                object.position.x = row[xChannel];
                object.position.y = row[yChannel];
                object.position.z = -1;
                setSelectElement(object);
                setCurrentAggregateSelectionFn(row);
              }
            }
          });
          if (!foundHex) {
            // if position is not inside a hexagon, we set select to null
            setSelectElement(null);
            setCurrentAggregateSelectionFn(null);
          }
        }
      } else {
        setSelectElement(null);
        setCurrentAggregateSelectionFn(null);
      }
      // eslint-disable-next-line
    }, [aggregateDataset, mouseClick, setCurrentAggregateSelectionFn]);

    return (
      <div>
        <LoadingIndicatorDialog
          handleClose={() => {
            cancelPromises();
          }}
          area={loadingArea}
        />
        {hexagons && hexagons.length > 0 && (
          <GLHexagons hexagons={hexagons} hoverElement={hoverElement} selectElement={selectElement} multipleId={multipleId} />
        )}
      </div>
    );
  },
);
