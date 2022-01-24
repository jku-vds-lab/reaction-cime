
import { useCancellablePromise } from 'projection-space-explorer'
import * as React from 'react'
import { connect, ConnectedProps } from "react-redux";
import * as THREE from 'three'
import { ReactionCIMEBackendFromEnv } from '../../Backend/ReactionCIMEBackend';
import { AppState } from '../../State/Store';
import { AggregateDataset } from '../AggregationTabPanel/AggregateDataset'
import { LoadingIndicatorDialog } from '../Dataset/DatasetTabPanel';
import { GLHeatmap } from './GLHeatmap'
import * as _ from "lodash";
import * as vsup from "vsup";
import * as d3 from 'd3v5';
import { setAggregateColorMapLegend } from '../../State/AggregateSettingsDuck';


const convert_to_rgb = (value: string | {r: number, g: number, b: number}):{r: number, g: number, b: number} => {
    if(Object.keys(value).includes("r") && Object.keys(value).includes("g") && Object.keys(value).includes("b"))
        return {"r": value["r"], "g": value["g"], "b": value["b"]};

    value = value.toString()
    if(value.startsWith("rgb")){
        var rgb = value.replace("rgb(", "")
        rgb = rgb.replace(")", "")
        rgb = rgb.replace(" ", "")
        var rgb_arr = rgb.split(",")
        return {"r": parseInt(rgb_arr[0]), "g": parseInt(rgb_arr[1]), "b": parseInt(rgb_arr[2])}
    }

    if(value.startsWith("#")){
        value = value.replace("#", "")
        var hex = value.match(/.{1,2}/g);
        return {"r": parseInt(hex[0], 16), "g": parseInt(hex[1], 16), "b": parseInt(hex[2], 16)}
    }
    
    console.log("error:", "format unknown ->", value)
    return {"r": 0, "g": 0, "b": 0};
}

const retrieve_information_from_agg_dataset = (aggregateDataset: AggregateDataset, dataset: AggregateDataset, value_col: string, uncertainty_col: string, setAggregateColorMapLegend: any, colorScale: String, useVSUP) => {
    if(!Object.keys(aggregateDataset.columns).includes(value_col)){
        return [null, null]
    }    
    var value_arr = dataset.vectors.map((row) => row[value_col]);

    // visit https://github.com/uwdata/vsup for more info about VSUP vs bivariate colorscale encoding
    var vDom = [aggregateDataset.columns[value_col].range.min, aggregateDataset.columns[value_col].range.max];

    var scale, legend, quantization;
    if(aggregateDataset.columns[uncertainty_col] != null){
        var uncertainty_arr = dataset.vectors.map((row) => row[uncertainty_col]);
        var uDom = [aggregateDataset.columns[uncertainty_col].range.min, aggregateDataset.columns[uncertainty_col].range.max];

        if(useVSUP){
            quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
            scale = vsup.scale().quantize(quantization).range(d3[colorScale]);
            
            legend = vsup.legend.arcmapLegend(scale);
            legend
                .vtitle(value_col)
                .utitle(uncertainty_col)
                // .size(180)
                // .x(20)
                // .y(50)
        }else{
            quantization = vsup.squareQuantization(4).valueDomain(vDom).uncertaintyDomain(uDom);

            scale = vsup.scale().quantize(quantization).range(d3[colorScale]);
            
            legend = vsup.legend.heatmapLegend(scale);
            legend
                .vtitle(value_col)
                .utitle(uncertainty_col)
                // .size(180)
                // .x(20)
                // .y(50)
        }
    }else{
        scale = d3.scaleQuantize()
            .domain(vDom)
            .range(d3.quantize(d3[colorScale], 8));

        legend = vsup.legend.simpleLegend()
            .title(value_col)
            // .height(20)
            .scale(scale)
            // .size(235)
            // .x(10)
            // .y(5)
    }

    // TODO: for extremely small values, we could add some scaling function that scales the values by some e-5 and add "e-5" to the label --> the ticks would then for example be "0.01" and the title "columnname e-5"
    // automatically obtain a suitable precision value
    const p = d3.precisionFixed((vDom[1]-vDom[0])/8);
    if(p > 2){
        legend.format(".2")
    }else if(p > 0){
        legend.format("." + p + "r");
    }

    setAggregateColorMapLegend(legend)

    // let background_colorMapping = new ContinuousMapping(
    //     {
    //         palette: [new SchemeColor('#fefefe'), new SchemeColor('#111111')],
    //         type: 'sequential'
    //     },
    //     aggregateDataset.columns[value_col].range // needs to have the same range as the overview Dataset
    // );



    var bgRGBA = new Uint8Array(value_arr.length * 4);
    for(var i=0; i < value_arr.length; i++){ // set opacity to 0 if value is not given
        if(value_arr[i] === undefined || isNaN(value_arr[i]) || value_arr[i] === ""){
            bgRGBA[4*i + 3] = 0;
        }else{

            var color;
            if(uncertainty_arr){
                color = scale(value_arr[i],uncertainty_arr[i])
            }else{
                color = scale(value_arr[i], 0)
            }
            var rgb = convert_to_rgb(color)
            
            // RGB from 0 to 255
            // var rgb = background_colorMapping.map(arr_pred[i]).rgb
            bgRGBA[4*i] = rgb.r;
            bgRGBA[4*i+1] = rgb.g;
            bgRGBA[4*i+2] = rgb.b;

            // OPACITY
            bgRGBA[4*i + 3] = 255;
        }
    }
    let bgDataTex = new THREE.DataTexture(bgRGBA, Math.sqrt(value_arr.length), Math.sqrt(value_arr.length), THREE.RGBAFormat);
    // bgDataTex.magFilter = THREE.LinearMipMapLinearFilter; // this causes border artefacts
    bgDataTex.magFilter = THREE.NearestFilter; // this makes it discrete
    // bgDataTex.magFilter = THREE.LinearFilter; // this causes border artefacts
    // bgDataTex.minFilter = THREE.LinearMipMapLinearFilter; // this causes everything to go black
    // bgDataTex.minFilter = THREE.LinearFilter;
    bgDataTex.minFilter = THREE.NearestFilter;


    let [width, height, x, y] = [100,100,0,0];
    if(dataset.bounds){
        width = dataset.bounds.x["max"]-dataset.bounds.x["min"];
        height = dataset.bounds.y["max"]-dataset.bounds.y["min"]; 

        // need to set x and y because the center is not at 0 0 if max and min are unequal
        x = (dataset.bounds.x["max"]+dataset.bounds.x["min"])/2;
        y = (dataset.bounds.y["max"]+dataset.bounds.y["min"])/2;
    }
    
    return [bgDataTex, {width:width, height:height, x:x, y:y}];

}

const mapStateToProps = (state: AppState) => ({
    aggregateColor: state.aggregateColor,
    poiDataset: state.dataset,
    viewTransform: state.viewTransform,
    colorScale: state.aggregateSettings?.colorscale,
    useVSUP: state.aggregateSettings?.useVSUP,
    sampleSize: state.aggregateSettings?.sampleSize
})
const mapDispatchToProps = (dispatch: any) => ({
    // setAggregateDataset: dataset => dispatch(setAggregateDatasetAction(dataset)),
    setAggregateColorMapLegend: (legend) => dispatch(setAggregateColorMapLegend(legend)),
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
}

const loading_area = "global_loading_indicator_aggregation_ds";
const AggregationLayer = connector(({ aggregateColor, poiDataset, viewTransform, setAggregateColorMapLegend, colorScale, useVSUP, sampleSize }: AggregationLayerProps) => {
    if(poiDataset == null || poiDataset.info == null || aggregateColor == null || aggregateColor.value_col == null || aggregateColor.value_col === "None"){
        return null;
    }

    let [textures, setTextures] = React.useState(null)
    let [sizes, setSizes] = React.useState(null)

    let [aggregateDataset, setAggregateDataset] = React.useState(null)
    let [aggregateDatasetZoomed, setAggregateDatasetZoomed] = React.useState(null)

    const { cancellablePromise, cancelPromises } = useCancellablePromise();

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

            // load zoomed version of the aggregate dataset
            ReactionCIMEBackendFromEnv.loadAggCSV((dataset) => {
                setAggregateDatasetZoomed(new AggregateDataset(dataset))
            }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateColor.cache_cols, sampleSize, range, cancellablePromise, abort_controller, "None")

        }, 500, {leading:false, trailing:true}) // leading: execute function at the beginning of the events; trailing: execute function at the end of the events; maxWait: maximum time the function is allowed to be delayed
        ,[aggregateColor, sampleSize]
    );

    React.useEffect(() => {
        if(aggregateDataset != null) // only start loading, if we already have a base dataset
            debouncedLoadAggDataset(viewTransform)
    // eslint-disable-next-line
    }, [viewTransform])

    React.useEffect(() => {
        // setAggregateDataset(null)
        let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
        // load the basic aggregateDataset with the high-level overview information
        ReactionCIMEBackendFromEnv.loadAggCSV((dataset) => {
            setAggregateDataset(new AggregateDataset(dataset))
        }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, aggregateColor.cache_cols, sampleSize, null, cancellablePromise, abort_controller, loading_area)

        // reset the zoomed version of the dataset
        setAggregateDatasetZoomed(null)
    // eslint-disable-next-line
    }, [aggregateColor, poiDataset.info.path, sampleSize])

    

    React.useEffect(() => {
        
        if(aggregateDataset && aggregateDataset.vectors){
            let sizes = [null, null];
            let textures = [null, null];
            let [texture, size] = retrieve_information_from_agg_dataset(aggregateDataset, aggregateDataset, aggregateColor.value_col, aggregateColor.uncertainty_col, setAggregateColorMapLegend, colorScale, useVSUP)
            sizes[0] = size;
            textures[0] = texture;

            if(aggregateDatasetZoomed && aggregateDatasetZoomed.vectors){
                let [texture, size] = retrieve_information_from_agg_dataset(aggregateDataset, aggregateDatasetZoomed, aggregateColor.value_col, aggregateColor.uncertainty_col, setAggregateColorMapLegend, colorScale, useVSUP)
                sizes[1] = size;
                textures[1] = texture;
            }
            setSizes(sizes);
            setTextures(textures);
        }else{
            setSizes(null);
            setTextures(null)
        }
    // eslint-disable-next-line
    }, [aggregateDataset, aggregateDatasetZoomed, colorScale, useVSUP])


    return <div>
    <LoadingIndicatorDialog
        handleClose={() => {
            cancelPromises();
        }}
        area={loading_area}
    />
    {(textures && textures.length > 0 && sizes && sizes.length > 0) && <GLHeatmap
            // texture={new THREE.TextureLoader().load('https://threejsfundamentals.org/threejs/resources/images/wall.jpg')}
            textures={textures}
            sizes={sizes}
        ></GLHeatmap>}
    </div>
    
})


export { AggregationLayer }