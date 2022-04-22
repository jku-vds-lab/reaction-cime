
import { useCancellablePromise } from 'projection-space-explorer'
import * as React from 'react'
import { connect, ConnectedProps } from "react-redux";
import * as THREE from 'three'
import { ReactionCIMEBackendFromEnv } from '../../../Backend/ReactionCIMEBackend';
import { AppState } from '../../../State/Store';
import { AggregateDataset } from '../AggregateDataset'
import { LoadingIndicatorDialog } from '../../Dataset/DatasetTabPanel';
import * as _ from "lodash";
import * as vsup from "vsup";
import * as d3 from 'd3v5';
import { GLContour } from './GLContour';
import { setAggregateColorMapScale } from '../../../State/AggregateSettingsDuck';

const retrieve_colorscale = (aggregateDataset: AggregateDataset, value_col: string, uncertainty_col: string, setAggregateColorMapScale: any, aggregateSettings: any) => {
    if(!Object.keys(aggregateDataset.columns).includes(value_col)){
        return [null, null]
    }    

    // visit https://github.com/uwdata/vsup for more info about VSUP vs bivariate colorscale encoding
    var vDom = [aggregateDataset.columns[value_col].range.min, aggregateDataset.columns[value_col].range.max];

    var scale, quantization;
    if(aggregateDataset.columns[uncertainty_col] != null){
        var uDom = [aggregateDataset.columns[uncertainty_col].range.min, aggregateDataset.columns[uncertainty_col].range.max];

        if(aggregateSettings?.useVSUP){
            quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
        }else{
            quantization = vsup.squareQuantization(4).valueDomain(vDom).uncertaintyDomain(uDom);
        }
        scale = vsup.scale().quantize(quantization).range(d3[aggregateSettings?.colorscale]);
    }else{
        scale = d3.scaleQuantize()
            .domain(vDom)
            .range(d3.quantize(d3[aggregateSettings.colorscale], 8));
    }
    setAggregateColorMapScale(scale)
    return scale;
}


const createContours = (dataset, value_col, scale) => {
    
    if(!Object.keys(dataset.columns).includes(value_col)){
        return [null, null]
    }
    // if(aggregateDataset.columns[uncertainty_col] != null){
    //     var uncertainty_arr = dataset.vectors.map((row) => row[uncertainty_col]);
    // }

    var value_arr = dataset.vectors.map((row) => row[value_col]);
    // console.log(value_arr)
    // let [width, height, x, y] = [100,100,0,0];
    // if(dataset.bounds){
    //     width = dataset.bounds.x["max"]-dataset.bounds.x["min"];
    //     height = dataset.bounds.y["max"]-dataset.bounds.y["min"]; 

    //     // need to set x and y because the center is not at 0 0 if max and min are unequal
    //     x = (dataset.bounds.x["max"]+dataset.bounds.x["min"])/2;
    //     y = (dataset.bounds.y["max"]+dataset.bounds.y["min"])/2;
    // }


    let xAxis = d3.scaleLinear()
        .range([0, Math.sqrt(value_arr.length)])
        .domain([dataset.bounds.x["min"], dataset.bounds.x["max"]])

    let yAxis = d3.scaleLinear()
        .range([0, Math.sqrt(value_arr.length)])
        .domain([dataset.bounds.y["min"], dataset.bounds.y["max"]])

    const min = Math.min(...value_arr)
    const max = Math.max(...value_arr)
    console.log(min, max)
    const noNanArr = value_arr.map((val)=> isNaN(parseFloat(val)) ? -Number.MAX_VALUE : val)
    console.log(noNanArr)
    let contours = d3.contours()
        // .thresholds(11)
        .thresholds(d3.range(min, max, (max-min)/11))
        .smooth(true)
        .size([Math.sqrt(value_arr.length), Math.sqrt(value_arr.length)])
        (noNanArr)
        

    let lines = []
    contours.forEach(contour => {
        // const coordinates = contour.coordinates[0][0]
        const coordinates_list = contour.coordinates
        const value = contour.value
        // let material = new THREE.LineBasicMaterial({ color: scale(value) })
        let material = new THREE.MeshBasicMaterial({ color: scale(value), side: THREE.DoubleSide })
        for (const key in coordinates_list) {
            const coordinates = coordinates_list[key][0]
            const points = []
            const shape = new THREE.Shape();
            shape.moveTo(xAxis.invert(coordinates[0][0]), yAxis.invert(coordinates[0][1]))
            for (let i = 0; i < coordinates.length ; i++) {
                let cur = coordinates[i]
                shape.lineTo(xAxis.invert(cur[0]), yAxis.invert(cur[1]))
                // points.push(new THREE.Vector2(xAxis.invert(cur[0]), yAxis.invert(cur[1])))

            }
            
            // var shape = new THREE.Shape()
            // shape.moveTo(xAxis.invert(coordinates[0][0]), yAxis.invert(coordinates[0][1]))
            // shape.splineThru(points)
            // var geo = new THREE.ShapeGeometry(shape);
            // let line = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(points), material)
            const geo = new THREE.ShapeGeometry(shape);
            let line = new THREE.Mesh(geo, material)


            line.visible = true
            lines.push(line)
        }

        
    })

    return lines
}

const mapStateToProps = (state: AppState) => ({
    aggregateColor: state.aggregateSettings?.aggregateColor,
    poiDataset: state.dataset,
    viewTransform: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.viewTransform,
    aggregateSettings: state.aggregateSettings,
})
const mapDispatchToProps = (dispatch: any) => ({
    // setAggregateDataset: dataset => dispatch(setAggregateDatasetAction(dataset)),
    setAggregateColorMapScale: (legend) => dispatch(setAggregateColorMapScale(legend)),
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
}

const loading_area = "global_loading_indicator_aggregation_ds";
export const AggregationContourLayer = connector(({ aggregateColor, poiDataset, viewTransform, setAggregateColorMapScale, aggregateSettings }: AggregationLayerProps) => {
    if(poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === "None"){
        return null;
    }

    let [lines, setLines] = React.useState(null)

    let [aggregateDataset, setAggregateDataset] = React.useState(null)

    const { cancellablePromise, cancelPromises } = useCancellablePromise();


    React.useEffect(() => {
        // setAggregateDataset(null)
        let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
        // load the basic aggregateDataset with the high-level overview information
        ReactionCIMEBackendFromEnv.loadAggCSV((dataset) => {
            setAggregateDataset(new AggregateDataset(dataset))
        }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateColor.cache_cols, aggregateSettings?.sampleSize, null, cancellablePromise, abort_controller, loading_area)

    }, [aggregateColor, poiDataset.info.path, aggregateSettings?.sampleSize])

    
    React.useEffect(() => {
        if(aggregateDataset && aggregateDataset.vectors){
            retrieve_colorscale(aggregateDataset, aggregateColor.value_col, aggregateColor.uncertainty_col, setAggregateColorMapScale, aggregateSettings) // "pred_var_9"
        }
    }, [aggregateDataset, aggregateColor, aggregateSettings?.colorscale, aggregateSettings?.useVSUP])

    React.useEffect(() => {
        if(aggregateSettings?.scale_obj != null && aggregateDataset && aggregateDataset.vectors){
            let lines = createContours(aggregateDataset, aggregateColor.value_col, aggregateSettings?.scale_obj) // "pred_var_9"
            setLines(lines);
        }else{
            setLines(null);
        }
    }, [aggregateSettings?.scale_obj, aggregateSettings?.valueFilter])


    return <div>
    <LoadingIndicatorDialog
        handleClose={() => {
            cancelPromises();
        }}
        area={loading_area}
    />
    {(lines && lines.length > 0) && <GLContour
            lines={lines}
        ></GLContour>}
    </div>
    
})


    
