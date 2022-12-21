import { RootState } from 'projection-space-explorer';
import React, { useState, useEffect } from 'react';
import * as THREE from 'three';
import { connect, ConnectedProps } from 'react-redux';

const mapStateToProps = (state: RootState) => ({
  viewTransform: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.viewTransform,
});

const mapDispatchToProps = (dispatch: any) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  /**
   * THREE texture
   */
  textures: [THREE.Texture];

  /**
   * World coordinate size
   */
  sizes: [
    {
      x: number;
      y: number;
      width: number;
      height: number;
    },
  ];
};

export const GLHeatmap = connector(({ viewTransform, textures, sizes }: Props) => {
  const ref = React.useRef<any>();

  const [renderer] = useState(
    () =>
      new THREE.WebGLRenderer({
        antialias: true,
        alpha: true,
      }),
  );
  const [scene] = useState(() => new THREE.Scene());
  const [rerender, setRerender] = useState(0);

  useEffect(() => {
    ref.current.appendChild(renderer.domElement);
    // eslint-disable-next-line
    }, [])

  useEffect(() => {
    scene.remove(...scene.children); // clears the scene
    for (const i in textures) {
      const texture = textures[i];
      const size = sizes[i];

      if (texture != null && size != null) {
        const groundMaterial = new THREE.MeshBasicMaterial({ map: texture, side: THREE.DoubleSide, transparent: true });
        const mesh = new THREE.Mesh(new THREE.PlaneBufferGeometry(size.width, size.height), groundMaterial);
        mesh.position.x = size.x;
        mesh.position.y = size.y;
        mesh.position.z = 0.0;
        // mesh.name="background" // could remove by name

        scene.add(mesh);
      }
    }
    setRerender(rerender + 1);
    // eslint-disable-next-line
    }, [textures, sizes])

  useEffect(() => {
    if (viewTransform) {
      const w = viewTransform.width;
      const h = viewTransform.height;

      renderer.setPixelRatio(window.devicePixelRatio);
      renderer.setSize(w, h);
      renderer.setClearColor(0x000000, 0);

      // Create orthographic camera
      const camera = new THREE.OrthographicCamera(w / -2, w / 2, h / 2, h / -2, 1, 1000);

      camera.position.z = 1;
      camera.lookAt(new THREE.Vector3(0, 0, 0));

      camera.position.x = viewTransform.centerX;
      camera.position.y = viewTransform.centerY;
      camera.zoom = viewTransform.zoom;

      camera.updateProjectionMatrix();

      try {
        renderer.render(scene, camera);
      } catch (e) {
        console.log(e);
      }
    }
    // eslint-disable-next-line
    }, [viewTransform, scene, scene.children, rerender])

  return <div style={{}} ref={ref} />;
});
