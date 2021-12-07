
import { CameraTransformations, ContinuousMapping, SchemeColor } from 'projection-space-explorer'
import * as React from 'react'
import { connect, ConnectedProps } from "react-redux";
import * as THREE from 'three'
import { AppState } from '../../State/Store';
import { AggregateDataset } from '../AggregationTabPanel/AggregateDataset'
import { GLHeatmap } from './GLHeatmap'



const mapStateToProps = (state: AppState) => ({
    aggregateDataset: state.aggregateDataset as AggregateDataset,
    aggregateColor: state.aggregateColor,
})
const mapDispatchToProps = (dispatch: any) => ({
    
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type AggregationLayerProps = PropsFromRedux & {
}


const AggregationLayer = connector(({ aggregateDataset }: AggregationLayerProps) => {
    
    let [texture, setTexture] = React.useState(null)
    let [width, setWidth] = React.useState(100);
    let [height, setHeight] = React.useState(100);

    React.useEffect(() => {

        if(aggregateDataset && aggregateDataset.vectors){
            
            setWidth(aggregateDataset.bounds.x["max"]-aggregateDataset.bounds.x["min"]);
            setHeight(aggregateDataset.bounds.y["max"]-aggregateDataset.bounds.y["min"]);

            var arr_pred = aggregateDataset.vectors.map((row) => row["val"]);
            // TODO: user input to choose colormap
            let background_colorMapping = new ContinuousMapping(
                {
                    palette: [new SchemeColor('#fefefe'), new SchemeColor('#111111')],
                    type: 'sequential'
                },
                // {palette:'dark2', type:'sequential'}, // TODO: use custom color scale
                // new ContinuosScale([
                //     new SchemeColor('#fefefe'),
                //     new SchemeColor('#111111')
                //   ]),
                aggregateDataset.columns["val"].range
            );

            var bgRGBA = new Uint8Array(arr_pred.length * 4);
            for(var i=0; i < arr_pred.length; i++){ // set opacity to 0 if value is not given
                if(arr_pred[i] === undefined || isNaN(arr_pred[i]) || arr_pred[i] === ""){
                    bgRGBA[4*i + 3] = 0;
                }else{
                    // RGB from 0 to 255
                    var rgb = background_colorMapping.map(arr_pred[i]).rgb
                    bgRGBA[4*i] = rgb.r;
                    bgRGBA[4*i+1] = rgb.g;
                    bgRGBA[4*i+2] = rgb.b;
                    // OPACITY
                    bgRGBA[4*i + 3] = 255;
                }
            }
            let bgDataTex = new THREE.DataTexture(bgRGBA, Math.sqrt(arr_pred.length), Math.sqrt(arr_pred.length), THREE.RGBAFormat);
            bgDataTex.magFilter = THREE.LinearMipMapLinearFilter;
            // bgDataTex.minFilter = THREE.LinearMipMapLinearFilter;
            
            setTexture(bgDataTex);
    
        }else{
            setTexture(null);
        }
    }, [aggregateDataset])
    
    return texture && <GLHeatmap
        // texture={new THREE.TextureLoader().load('https://threejsfundamentals.org/threejs/resources/images/wall.jpg')}
        texture={texture}
        size={{
        x: 0,
        y: 0,
        width: width,
        height: height
        }}
    ></GLHeatmap>
    
})


export { AggregationLayer }