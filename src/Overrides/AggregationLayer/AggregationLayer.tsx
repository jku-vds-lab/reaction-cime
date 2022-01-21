
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



const mapStateToProps = (state: AppState) => ({
    aggregateColor: state.aggregateColor,
    poiDataset: state.dataset,
    viewTransform: state.viewTransform
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
const AggregationLayer = connector(({ aggregateColor, poiDataset, viewTransform, setAggregateColorMapLegend }: AggregationLayerProps) => {
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
            }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, range, cancellablePromise, abort_controller, "None")

        }, 500, {leading:false, trailing:true}) // leading: execute function at the beginning of the events; trailing: execute function at the end of the events; maxWait: maximum time the function is allowed to be delayed
        ,[aggregateColor]
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
        }, poiDataset.info.path, aggregateColor.value_col, aggregateColor.uncertainty_col, null, cancellablePromise, abort_controller, loading_area)

        // reset the zoomed version of the dataset
        setAggregateDatasetZoomed(null)
    // eslint-disable-next-line
    }, [aggregateColor, poiDataset.info.path])

    const retrieve_information_from_agg_dataset = (dataset: AggregateDataset) => {
        var value_arr = dataset.vectors.map((row) => row["val"]);
        var vDom = [aggregateDataset.columns["val"].range.min, aggregateDataset.columns["val"].range.max];

        var scale, legend;
        if(aggregateDataset.columns["uncertainty"] != null){
            var uncertainty_arr = dataset.vectors.map((row) => row["uncertainty"]);
            var uDom = [aggregateDataset.columns["uncertainty"].range.min, aggregateDataset.columns["uncertainty"].range.max];
            
            var quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
            
            scale = vsup.scale().quantize(quantization).range(d3.interpolateYlGnBu);
            
            legend = vsup.legend.arcmapLegend(scale);
            legend
                .vtitle(aggregateColor.value_col)
                .utitle(aggregateColor.uncertainty_col)
                .size(180)
                .x(20)
                .y(50)
        }else{
            scale = d3.scaleQuantize()
                .domain(vDom)
                .range(d3.quantize(d3.interpolateYlGnBu, 8));

            legend = vsup.legend.simpleLegend()
                .title(aggregateColor.value_col)
                .height(20)
                .scale(scale)
                .size(235)
                .x(10)
                .y(5)
        }

        // TODO: for extremely small values, we could add some scaling function that scales the values by some e-5 and add "e-5" to the label --> the ticks would then for example be "0.01" and the title "columnname e-5"
        // automatically obtain a suitable precision value
        const p = d3.precisionFixed((vDom[1]-vDom[0])/8);
        if(p > 2){
            legend.format(".2")
        }else{
            legend.format("." + p + "r");
        }

        setAggregateColorMapLegend(legend)

        // let background_colorMapping = new ContinuousMapping(
        //     {
        //         palette: [new SchemeColor('#fefefe'), new SchemeColor('#111111')],
        //         type: 'sequential'
        //     },
        //     aggregateDataset.columns["val"].range // needs to have the same range as the overview Dataset
        // );



        var bgRGBA = new Uint8Array(value_arr.length * 4);
        for(var i=0; i < value_arr.length; i++){ // set opacity to 0 if value is not given
            if(value_arr[i] === undefined || isNaN(value_arr[i]) || value_arr[i] === ""){
                bgRGBA[4*i + 3] = 0;
            }else{
                // RGB from 0 to 255
                // var rgb = background_colorMapping.map(arr_pred[i]).rgb
                // bgRGBA[4*i] = rgb.r;
                // bgRGBA[4*i+1] = rgb.g;
                // bgRGBA[4*i+2] = rgb.b;

                var rgb;
                if(uncertainty_arr){
                    rgb = scale(value_arr[i],uncertainty_arr[i])
                }else{
                    rgb = scale(value_arr[i], 0)
                }
                rgb = rgb.replace("rgb(", "")
                rgb = rgb.replace(")", "")
                rgb = rgb.replace(" ", "")
                rgb = rgb.split(",")
                bgRGBA[4*i] = rgb[0];
                bgRGBA[4*i+1] = rgb[1];
                bgRGBA[4*i+2] = rgb[2];

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

    React.useEffect(() => {
        
        if(aggregateDataset && aggregateDataset.vectors){
            let sizes = [null, null];
            let textures = [null, null];
            let [texture, size] = retrieve_information_from_agg_dataset(aggregateDataset)
            sizes[0] = size;
            textures[0] = texture;

            if(aggregateDatasetZoomed && aggregateDatasetZoomed.vectors){
                let [texture, size] = retrieve_information_from_agg_dataset(aggregateDatasetZoomed)
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
    }, [aggregateDataset, aggregateDatasetZoomed])


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