import { RootState } from "projection-space-explorer";
import React, { useState } from "react";
import { useEffect } from "react";
import * as THREE from 'three';
import { connect, ConnectedProps } from "react-redux";


const mapStateToProps = (state: RootState) => ({
    viewTransform: state.viewTransform
})

const mapDispatchToProps = (dispatch: any) => ({
})


const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>

type Props = PropsFromRedux & {
    /**
     * THREE texture
     */
    lines: [THREE.Mesh]
    
}


export const GLContour = connector(({ viewTransform, lines }: Props) => {
    const ref = React.useRef<any>()

    const [renderer] = useState(() => new THREE.WebGLRenderer({
        antialias: true,
        alpha: true
    }))
    const [scene] = useState(() => new THREE.Scene())
    const [rerender, setRerender] = useState(0)

    useEffect(() =>  {
        ref.current.appendChild(renderer.domElement);
        // eslint-disable-next-line
    }, [])

    useEffect(() => {
        scene.clear(); // clears the scene //.remove(...scene.children)//
        for(let i in lines){
            const line = lines[i]
            
            if(lines != null){
                scene.add(line);
            }
        }
        setRerender(rerender+1)
        // eslint-disable-next-line
    }, [lines])


    useEffect(() => {
        if (viewTransform) {
            const w = viewTransform.width
            const h = viewTransform.height

            renderer.setPixelRatio(window.devicePixelRatio);
            renderer.setSize(w, h)
            renderer.setClearColor( 0x000000, 0 );

            // Create orthographic camera
            var camera = new THREE.OrthographicCamera(w / - 2, w / 2, h / 2, h / - 2, 1, 1000);

            camera.position.z = 1;
            camera.lookAt(new THREE.Vector3(0, 0, 0));

            camera.position.x = viewTransform.centerX
            camera.position.y = viewTransform.centerY
            camera.zoom = viewTransform.zoom

            camera.updateProjectionMatrix();

            try {
                renderer.render(scene, camera);
            } catch (e) {
                console.log(e)
            }
        }
    // eslint-disable-next-line
    }, [viewTransform, scene, scene.children, rerender])

    return <div style={{  }} ref={ref} />
})