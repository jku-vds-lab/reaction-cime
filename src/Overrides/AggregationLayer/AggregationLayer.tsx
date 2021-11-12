
import { CameraTransformations, RenderingContextEx } from 'projection-space-explorer'
import * as React from 'react'
import { connect } from 'react-redux'
import * as THREE from 'three'
import { GLHeatmap } from './GLHeatmap'

type AggregationLayerProps = {
    viewTransform: CameraTransformations
}

const mapStateToProps = state => ({
    viewTransform: state.viewTransform as CameraTransformations
})
const mapDispatchToProps = (dispatch: any) => ({
})

const connector = connect(mapStateToProps, mapDispatchToProps);

const AggregationLayer = connector(({  }: AggregationLayerProps) => {
    var arr = [[-68.63237573405219,-72.60442862855822,]
        ,[-33.4210922694048,-72.60442862855822,]
        ,[1.7901911952425849,-72.60442862855822,0.0]
        ,[37.00147465988998,-72.60442862855822,]
        ,[72.21275812453736,-72.60442862855822,]
        ,[-68.63237573405219,-35.970975946614715,]
        ,[-33.4210922694048,-35.970975946614715,9.230136385552946]
        ,[1.7901911952425849,-35.970975946614715,16.684928472095358]
        ,[37.00147465988998,-35.970975946614715,22.972915986783004]
        ,[72.21275812453736,-35.970975946614715,]
        ,[-68.63237573405219,0.662476735328795,]
        ,[-33.4210922694048,0.662476735328795,22.018889439919846]
        ,[1.7901911952425849,0.662476735328795,17.09708391241856]
        ,[37.00147465988998,0.662476735328795,21.170857542513634]
        ,[72.21275812453736,0.662476735328795,]
        ,[-68.63237573405219,37.295929417272305,]
        ,[-33.4210922694048,37.295929417272305,9.370991882020375]
        ,[1.7901911952425849,37.295929417272305,9.561317447766088]
        ,[37.00147465988998,37.295929417272305,6.001010531532636]
        ,[72.21275812453736,37.295929417272305,]
        ,[-68.63237573405219,73.9293820992158,]
        ,[-33.4210922694048,73.9293820992158,]
        ,[1.7901911952425849,73.9293820992158,]
        ,[37.00147465988998,73.9293820992158,]
        ,[72.21275812453736,73.9293820992158,]];


    var arr_pred = arr.map((row) => row[2]);

    var dummyRGBA = new Uint8Array(arr_pred.length * 4);
    for(var i=0; i < arr_pred.length; i++){
        console.log(arr_pred[i])
        // RGB from 0 to 255
        dummyRGBA[4*i] = 255*arr_pred[i]/23;
        // OPACITY
        dummyRGBA[4*i + 3] = 255;
        // if(arr_pred[i] === undefined){
        //     dummyRGBA[4*i + 3] = 0;
        // }
    }
    let dummyDataTex = new THREE.DataTexture(dummyRGBA, Math.sqrt(arr_pred.length), Math.sqrt(arr_pred.length), THREE.RGBAFormat);
    dummyDataTex.magFilter = THREE.LinearMipMapLinearFilter;
    // dummyDataTex.minFilter = THREE.LinearMipMapLinearFilter;

    let texture = dummyDataTex//new THREE.TextureLoader().load('https://threejsfundamentals.org/threejs/resources/images/wall.jpg');
    console.log(texture)
    return <GLHeatmap
    texture={texture}
    size={{
      x: 0,
      y: 0,
      width: 100,
      height: 100
    }}
  ></GLHeatmap>
    
})


export { AggregationLayer }