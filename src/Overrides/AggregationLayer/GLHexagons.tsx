import { RootState } from 'projection-space-explorer';
import React, { useState, useEffect } from 'react';
import * as THREE from 'three';
import { connect, ConnectedProps } from 'react-redux';
import { EntityId } from '@reduxjs/toolkit';

const mapStateToProps = (state: RootState) => ({
  smallMultiples: state.multiples.multiples.entities,
});

const mapDispatchToProps = (dispatch: any) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  /**
   * THREE texture
   */
  hexagons: [THREE.Mesh];
  hoverElement: THREE.Mesh;
  selectElement: THREE.Mesh;
  multipleId: EntityId;
};

export const GLHexagons = connector(({ smallMultiples, hexagons, hoverElement, selectElement, multipleId }: Props) => {
  const { viewTransform } = smallMultiples[multipleId].attributes;
  const ref = React.useRef<any>();
  const canvasRef = React.useRef<HTMLCanvasElement>();

  const [dim, setDim] = useState({ width: 0, height: 0 });

  const [renderer, setRenderer] = useState<THREE.WebGLRenderer>();
  const [scene] = useState(() => new THREE.Scene());
  const [rerender, setRerender] = useState(0);

  useEffect(() => {
    // ref.current.appendChild(renderer.domElement);
    // eslint-disable-next-line
    setRenderer(
      new THREE.WebGLRenderer({
        antialias: true,
        alpha: true,
        canvas: canvasRef.current,
      }),
    );
  }, []);

  useEffect(() => {
    scene.remove(...scene.children);
    hexagons?.forEach((hex) => {
      scene.add(hex);
    });
    setRerender(rerender + 1);
    // eslint-disable-next-line
  }, [hexagons]);

  useEffect(() => {
    if (!renderer) return () => {};

    const container = canvasRef.current;

    const observer = new ResizeObserver(() => {
      const w = container.offsetWidth;
      const h = container.offsetHeight;

      setDim({ width: w, height: h });

      renderer?.setSize(container.offsetWidth, container.offsetHeight, false);

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
    });

    observer.observe(container);

    return () => {
      observer.unobserve(container);
    };
  }, [renderer, canvasRef, scene, viewTransform]);

  useEffect(() => {
    const selectedObject = scene.getObjectByName('hoverElement');
    scene.remove(selectedObject);

    if (hoverElement != null) {
      hoverElement.name = 'hoverElement';
      scene.add(hoverElement);
    }
    setRerender(rerender + 1);
    // eslint-disable-next-line
  }, [hoverElement]);

  useEffect(() => {
    const selectedObject = scene.getObjectByName('selectElement');
    scene.remove(selectedObject);

    if (selectElement != null) {
      selectElement.name = 'selectElement';
      scene.add(selectElement);
    }
    setRerender(rerender + 1);
    // eslint-disable-next-line
  }, [selectElement]);

  useEffect(() => {
    if (viewTransform && renderer) {
      const w = viewTransform.width;
      const h = viewTransform.height;

      renderer.setPixelRatio(window.devicePixelRatio);
      // renderer.setSize(w - 16, h - 16);
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
  }, [viewTransform, scene, scene.children, rerender]);

  return <canvas ref={canvasRef} style={{ width: '100%', height: '100%' }} width={dim.width} height={dim.height} />;
});
