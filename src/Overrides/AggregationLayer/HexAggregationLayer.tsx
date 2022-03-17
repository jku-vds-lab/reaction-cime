
import { useCancellablePromise } from 'projection-space-explorer'
import * as React from 'react'
import { connect, ConnectedProps } from "react-redux";
import * as THREE from 'three'
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { AppState } from '../../State/Store';
import { AggregateDataset } from './AggregateDataset'
import { LoadingIndicatorDialog } from '../Dataset/DatasetTabPanel';
import { setUncertaintyRange, setValueRange } from '../../State/AggregateSettingsDuck';
import { GLHexagons } from './GLHexagons';
import * as _ from "lodash";
import { PSE_BLUE } from '../../Utility/Utils';
import { setCurrentAggregateSelection } from '../../State/SelectionDuck';


const createHexagons = (dataset: AggregateDataset, value_col: string, uncertainty_col: string, scale_obj:any, valueFilter: string[]) => {
    
    let hexagons = []

    dataset.vectors.forEach((row) => {
        if(row["hex"] === "False"){
            console.log("wrong point")
            let materialWrong = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.DoubleSide })
            var geometryWrong = new THREE.CircleGeometry(1, 100)
            var objectWrong = new THREE.Mesh(geometryWrong, materialWrong)
            objectWrong.position.x = row.x
            objectWrong.position.y = row.y

            hexagons.push(objectWrong)
        }else{
            if(row[value_col] === undefined || isNaN(row[value_col]) || row[value_col] === ""){

            }else{
                var color;
                if(dataset.columns[uncertainty_col] != null){
                    color = scale_obj(row[value_col],row[uncertainty_col])
                }else{
                    color = scale_obj(row[value_col], 0)
                }
                if(valueFilter == null || valueFilter.length <= 0 || valueFilter.includes(color)){ // if there is a filter, and the current value is not within, we ignore this hexagon
                    let material = new THREE.MeshBasicMaterial({ color: color, side: THREE.DoubleSide })
                    var geometry = new THREE.CircleGeometry(row["circ_radius"]-row["circ_radius"]*0.04, 6)
                    var object = new THREE.Mesh(geometry, material)
                    object.position.x = row.x
                    object.position.y = row.y
        
                    hexagons.push(object)
                }
            }
        }
        
    })

    return hexagons
}

// test if point is inside hexagon adapted from http://www.playchilla.com/how-to-check-if-a-point-is-inside-a-hexagon
function isInsideHex(point_x: number, point_y: number, hex_x: number, hex_y: number, radius: number, circ_radius: number){
    const point_q2x = Math.abs(point_x-hex_x) // transform the test point locally and to quadrant 2
    const point_q2y = Math.abs(point_y-hex_y) // transform the test point locally and to quadrant 2
    if((point_q2x > circ_radius) || (point_q2y > radius)){
        return false;
    }
    return (radius*circ_radius - radius*point_q2x - circ_radius/2*point_q2y) >= 0 // finally the dot product can be reduced to this due to the hexagon symmetry
    
}

// test if point is outside of a bounding box, including a certain threshold
function isOutsideBoundingBox(position, dataset, threshold){
    return position.x < (dataset.columns["x"].range.min - threshold) 
                || position.x > (dataset.columns["x"].range.max + threshold)
                || position.y < (dataset.columns["y"].range.min - threshold) 
                || position.y > (dataset.columns["y"].range.max + threshold)
}
    

const mapStateToProps = (state: AppState) => ({
    aggregateColor: state.aggregateSettings?.aggregateColor,
    poiDataset: state.dataset,
    viewTransform: state.viewTransform,
    aggregateSettings: state.aggregateSettings,
    mouseMove: state.mouseInteractionHooks?.mousemove,
    mouseClick: state.mouseInteractionHooks?.mouseclick,
})
const mapDispatchToProps = (dispatch: any) => ({
    setValueRange: (range) => dispatch(setValueRange(range)),
    setUncertaintyRange: (range) => dispatch(setUncertaintyRange(range)),
    setCurrentAggregateSelectionFn: (selection) => dispatch(setCurrentAggregateSelection(selection))
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
}

const loading_area = "global_loading_indicator_aggregation_ds";
export const HexAggregationLayer = connector(({ setCurrentAggregateSelectionFn, aggregateColor, poiDataset, viewTransform, setValueRange, setUncertaintyRange, aggregateSettings, mouseMove, mouseClick }: AggregationLayerProps) => {
    
    if(poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === "None"){
        return null;
    }

    const [hexagons, setHexagons] = React.useState(null)
    const [hoverElement, setHoverElement] = React.useState(null)
    const [selectElement, setSelectElement] = React.useState(null)
    const [aggregateDataset, setAggregateDataset] = React.useState<AggregateDataset>(null)
    const [datasetValueRange, setDatasetValueRange] = React.useState(null)
    const [datasetUncertaintyRange, setDatasetUncertaintyRange] = React.useState(null)

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    React.useEffect(() => {
        if(aggregateColor?.value_col != null){
            ReactionCIMEBackendFromEnv.loadValueRange(poiDataset.info.path, aggregateColor.value_col).then((response) => {
                setDatasetValueRange(response)
            })
            if(aggregateColor?.uncertainty_col != null){
                ReactionCIMEBackendFromEnv.loadValueRange(poiDataset.info.path, aggregateColor.uncertainty_col).then((response) => {
                    setDatasetUncertaintyRange(response)
                })
            }else{
                setDatasetUncertaintyRange(null)
            }
        }
    }, [aggregateColor, poiDataset.info.path])
    

    React.useEffect(() => {
        if(aggregateSettings?.deriveRange){
            if(datasetValueRange != null){
                setValueRange(datasetValueRange)
                setUncertaintyRange(datasetUncertaintyRange)
            }
        }
        // eslint-disable-next-line
    }, [datasetValueRange, datasetUncertaintyRange, aggregateSettings?.deriveRange])
    
    // eslint-disable-next-line react-hooks/exhaustive-deps
const debouncedLoadAggDataset = React.useCallback(
    _.debounce((viewTransform) => {
        cancelPromises();
        let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
        
        const range = {
            x: {
                min: Math.round(-viewTransform.width/viewTransform.zoom/2 + viewTransform.centerX), 
                max: Math.round(viewTransform.width/viewTransform.zoom/2 + viewTransform.centerX)
            },
            y: {
                min: Math.round(-viewTransform.height/viewTransform.zoom/2 + viewTransform.centerY), 
                max: Math.round(viewTransform.height/viewTransform.zoom/2 + viewTransform.centerY)
            }
        }
        
        // take 150% of the boundaries, such that we have clean borders
        range.x.min = range.x.min - Math.abs(range.x.min/2)
        range.x.max = range.x.max + Math.abs(range.x.max/2)
        range.y.min = range.y.min - Math.abs(range.y.min/2)
        range.y.max = range.y.max + Math.abs(range.y.max/2)

        // load the basic aggregateDataset with the high-level overview information
        ReactionCIMEBackendFromEnv.loadHexAgg((dataset) => {
            setAggregateDataset(new AggregateDataset(dataset))
        }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateColor.cache_cols, aggregateSettings?.sampleSize, aggregateSettings?.aggregationMethod, range, cancellablePromise, abort_controller, loading_area)
        
    }, 500, {leading:false, trailing:true}) // leading: execute function at the beginning of the events; trailing: execute function at the end of the events; maxWait: maximum time the function is allowed to be delayed
    ,[poiDataset?.info?.path, aggregateColor, aggregateSettings]
);

    React.useEffect(() => {
        // setAggregateDataset(null)
        debouncedLoadAggDataset(viewTransform)
    // eslint-disable-next-line
    }, [aggregateColor, poiDataset.info.path, aggregateSettings?.sampleSize, aggregateSettings?.aggregationMethod, viewTransform])
    

    React.useEffect(() => {
        if(aggregateSettings?.scale_obj != null){
            if(aggregateDataset && aggregateDataset.vectors){
                if(Object.keys(aggregateDataset.columns).includes(aggregateColor.value_col)){
                    
                    let hexs = createHexagons(aggregateDataset, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateSettings?.scale_obj, aggregateSettings?.valueFilter) // "pred_var_9"
                    setHexagons(hexs);
                }
            }else{
                setHexagons(null);
            }
        }
        
        // aggregateColor --> aggregateColor has direct influence on aggregateDataset through "setAggregateDataset" in the useEffect above
        // eslint-disable-next-line
    }, [aggregateSettings?.scale_obj, aggregateDataset, aggregateSettings?.valueFilter])


    React.useEffect(() => {
        if(mouseMove != null && !mouseMove.event_used && aggregateDataset != null){
            const circ_radius = aggregateDataset.columns["circ_radius"].range.max
            // check if it is outside of the aggregation bounding box
            if(isOutsideBoundingBox(mouseMove, aggregateDataset, circ_radius)){
                setHoverElement(null)
            }else{ // if it is inside, we need to check each hexagon, if coordinates are inside
                const radius = Math.sqrt(3)/2 * circ_radius;// https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
                let foundHex = false;
                aggregateDataset.vectors.forEach((row) => {
                    const isInside = isInsideHex(mouseMove.x, mouseMove.y, row.x, row.y, radius, circ_radius)
                    if(isInside){
                        foundHex = true;
                        if(hoverElement == null || hoverElement.position.x !== row.x || hoverElement.position.y !== row.y){
                            let material = new THREE.MeshBasicMaterial({ color: PSE_BLUE, side: THREE.DoubleSide, transparent: true, opacity: 1.0 })
                            var geometry = new THREE.CircleGeometry(row["circ_radius"]-row["circ_radius"]*0.04, 6)
                            var object = new THREE.Mesh(geometry, material)
                            object.position.x = row.x
                            object.position.y = row.y
                            setHoverElement(object)
                        }
                        return;
                    }
                })

                if(!foundHex){ // if position is not inside a hexagon, we set hover to null
                    setHoverElement(null)
                }
            }
        }else{
            setHoverElement(null)
        }
        
        // eslint-disable-next-line
    }, [aggregateDataset, mouseMove])

    React.useEffect(() => {
        if(mouseClick != null && !mouseClick.event_used && aggregateDataset != null){
            const circ_radius = aggregateDataset.columns["circ_radius"].range.max
            // check if it is outside of the aggregation bounding box
            if(isOutsideBoundingBox(mouseClick, aggregateDataset, circ_radius)){
                setSelectElement(null)
                setCurrentAggregateSelectionFn(null)
            }else{ // if it is inside, we need to check each hexagon, if coordinates are inside
                const radius = Math.sqrt(3)/2 * circ_radius;// https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
                let foundHex = false;
                aggregateDataset.vectors.forEach((row) => {
                    const isInside = isInsideHex(mouseClick.x, mouseClick.y, row.x, row.y, radius, circ_radius)
                    if(isInside){
                        foundHex = true;
                        if(selectElement == null || selectElement.position.x !== row.x || selectElement.position.y !== row.y){
                            let material = new THREE.MeshBasicMaterial({ color: PSE_BLUE, side: THREE.DoubleSide, transparent: true, opacity: 1.0 })
                            const radius = Number(row["circ_radius"])+Number(row["circ_radius"])*0.04
                            var geometry = new THREE.CircleGeometry(radius, 6)
                            var object = new THREE.Mesh(geometry, material)
                            object.position.x = row.x
                            object.position.y = row.y
                            object.position.z = -1
                            setSelectElement(object)
                            console.log(row)
                            setCurrentAggregateSelectionFn(row)
                        }
                        return;
                    }
                })
                if(!foundHex){ // if position is not inside a hexagon, we set select to null
                    setSelectElement(null)
                    setCurrentAggregateSelectionFn(null)
                }
            }
        }
        
    }, [aggregateDataset, mouseClick, setCurrentAggregateSelectionFn])


    return <div>
    <LoadingIndicatorDialog
        handleClose={() => {
            cancelPromises();
        }}
        area={loading_area}
    />
    {(hexagons && hexagons.length > 0) && <GLHexagons
            hexagons={hexagons}
            hoverElement={hoverElement}
            selectElement={selectElement}
        ></GLHexagons>}
    </div>
    
})


    
