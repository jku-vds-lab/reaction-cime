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
    texture: THREE.Texture

    /**
     * World coordinate size
     */
    size: {
        x: number
        y: number
        width: number
        height: number
    }
}


export const GLHeatmap = connector(({ viewTransform, texture, size }: Props) => {
    const ref = React.useRef<any>()

    const [renderer] = useState(() => new THREE.WebGLRenderer({
        antialias: true,
        alpha: true
    }))
    const [scene] = useState(() => new THREE.Scene())
    const [mesh, setMesh] = useState<THREE.Mesh>(new THREE.Mesh())


    useEffect(() => {
        ref.current.appendChild(renderer.domElement);

        // var geometry = new THREE.BoxGeometry(16, 16, 16);
        // var material = new THREE.MeshBasicMaterial({ color: 0x00ff00 });

        var groundMaterial = new THREE.MeshBasicMaterial({ map: texture, side: THREE.DoubleSide });
        var mesh = new THREE.Mesh(new THREE.PlaneBufferGeometry(size.width, size.height), groundMaterial);
        mesh.position.y = size.x;
        mesh.position.x = size.y;
        mesh.position.z = 0.0;

        // const cube = new THREE.Mesh(geometry, material);
        scene.add(mesh);

        setMesh(mesh)
    }, [])

    useEffect(() => {

    }, [texture])

    useEffect(() => {
        if (viewTransform) {
            const w = viewTransform.width
            const h = viewTransform.height

            renderer.setPixelRatio(window.devicePixelRatio);
            renderer.setSize(w, h)

            // Create orthographic camera
            var camera = new THREE.OrthographicCamera(w / - 2, w / 2, h / 2, h / - 2, 1, 1000);

            camera.position.z = 100;
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
    }, [viewTransform, mesh])

    return <div style={{  }} ref={ref} />
})