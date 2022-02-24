
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


const createHexagons = (dataset: AggregateDataset, value_col: string, uncertainty_col: string, scale_obj:any, valueFilter: string[]) => {
    
    let hexagons = []

    dataset.vectors.forEach((row) => {
        if(row["hex"] === "False"){
            console.log("wrong point")
            let material = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.DoubleSide })
            var geometry = new THREE.CircleGeometry(1, 100)
            var object = new THREE.Mesh(geometry, material)
            object.position.x = row.x
            object.position.y = row.y

            hexagons.push(object)
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
    

const mapStateToProps = (state: AppState) => ({
    aggregateColor: state.aggregateSettings?.aggregateColor,
    poiDataset: state.dataset,
    viewTransform: state.viewTransform,
    aggregateSettings: state.aggregateSettings,
    mouseMove: state.mouseInteractionHooks?.mousemove,
})
const mapDispatchToProps = (dispatch: any) => ({
    setValueRange: (range) => dispatch(setValueRange(range)),
    setUncertaintyRange: (range) => dispatch(setUncertaintyRange(range)),
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
}

const loading_area = "global_loading_indicator_aggregation_ds";
export const HexAggregationLayer = connector(({ aggregateColor, poiDataset, viewTransform, setValueRange, setUncertaintyRange, aggregateSettings, mouseMove }: AggregationLayerProps) => {
    if(poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === "None"){
        return null;
    }

    let [hexagons, setHexagons] = React.useState(null)
    let [hoverElement, setHoverElement] = React.useState(null)
    let [aggregateDataset, setAggregateDataset] = React.useState<AggregateDataset>(null)

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

    
    React.useEffect(() => {
        // setAggregateDataset(null)
        let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
        // load the basic aggregateDataset with the high-level overview information
        ReactionCIMEBackendFromEnv.loadHexAgg((dataset) => {
            setAggregateDataset(new AggregateDataset(dataset))
        }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateColor.cache_cols, aggregateSettings?.sampleSize, aggregateSettings?.aggregationMethod, cancellablePromise, abort_controller, loading_area)

    }, [aggregateColor, poiDataset.info.path, aggregateSettings?.sampleSize, aggregateSettings?.aggregationMethod])


    React.useEffect(() => {
        if(aggregateDataset && aggregateDataset.vectors && aggregateSettings?.deriveRange){
            if(Object.keys(aggregateDataset.columns).includes(aggregateColor.value_col)){
                setValueRange(aggregateDataset.columns[aggregateColor.value_col].range)
            
                if(aggregateDataset.columns[aggregateColor.uncertainty_col] != null){
                    setUncertaintyRange(aggregateDataset.columns[aggregateColor.uncertainty_col].range)
                }else{
                    setUncertaintyRange(null)
                }
            }
        }
        // aggregateColor --> aggregateColor has direct influence on aggregateDataset through "setAggregateDataset" in the useEffect above
        // eslint-disable-next-line
    }, [aggregateDataset, aggregateSettings?.deriveRange])


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
        
    }, [aggregateSettings?.scale_obj, aggregateDataset, aggregateSettings?.valueFilter])


    React.useEffect(() => {
        if(mouseMove != null && !mouseMove.event_used && aggregateDataset != null){
            const circ_radius = aggregateDataset.columns["circ_radius"].range.max
            // check if it is outside of the aggregation bounding box
            if(mouseMove.x < (aggregateDataset.columns["x"].range.min - circ_radius) 
                || mouseMove.x > (aggregateDataset.columns["x"].range.max + circ_radius)
                || mouseMove.y < (aggregateDataset.columns["y"].range.min - circ_radius) 
                || mouseMove.y > (aggregateDataset.columns["y"].range.max + circ_radius)){
                setHoverElement(null)
            }else{ // if it is inside, we need to check each hexagon, if coordinates are inside
                const radius = Math.sqrt(3)/2 * circ_radius;// https://www-formula.com/geometry/circle-inscribed/radius-circle-inscribed-regular-hexagon
                let foundHex = false;
                aggregateDataset.vectors.forEach((row) => {
                    const isInside = isInsideHex(mouseMove.x, mouseMove.y, row.x, row.y, radius, circ_radius)
                    if(isInside){
                        foundHex = true;
                        if(hoverElement == null || hoverElement.position.x !== row.x || hoverElement.position.y !== row.y){
                            let material = new THREE.MeshBasicMaterial({ color: 0x000000, side: THREE.DoubleSide, transparent: true, opacity: 0.5 })
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
    }, [aggregateDataset, mouseMove])

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
        ></GLHexagons>}
    </div>
    
})


    
