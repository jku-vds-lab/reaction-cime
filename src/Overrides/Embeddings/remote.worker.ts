/* eslint-disable no-restricted-globals */
// import "regenerator-runtime/runtime";

import { RemoteEmbedding } from './RemoteEmbedding';

const callbackFnStep = (emb, step, msg) => {
  self.postMessage({ messageType: 'step', embedding: emb, step, msg });
};
const callbackFnTerminated = () => {
  self.postMessage({ messageType: 'terminated' });
};

let embedding: RemoteEmbedding | null;

self.addEventListener(
  'message',
  (e) => {
    if (e.data.messageType === 'init') {
      embedding = new RemoteEmbedding(e.data.backendUrl, e.data.path, e.data.params, e.data.init_coordinates, e.data.selected_feature_info, e.data.ids);
      embedding.initializeFit(callbackFnStep);
    } else if (e.data.messageType === 'abort') {
      embedding?.abort(callbackFnTerminated);
    }
  },
  false,
);

const helper = null;
export default helper;
