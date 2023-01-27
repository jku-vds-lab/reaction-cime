// import "regenerator-runtime/runtime";

import { RemoteEmbedding } from './RemoteEmbedding';

/* eslint-disable-next-line no-restricted-globals */
const ctx = self;

const callbackFnStep = (emb, step, msg) => {
  ctx.postMessage({ messageType: 'step', embedding: emb, step, msg });
};
const callbackFnTerminated = () => {
  ctx.postMessage({ messageType: 'terminated' });
};

ctx.addEventListener(
  'message',
  function (e) {
    if (e.data.messageType === 'init') {
      const embedding = new RemoteEmbedding(e.data.backendUrl, e.data.path, e.data.params, e.data.init_coordinates, e.data.selected_feature_info);
      embedding.initializeFit(callbackFnStep);
      // embedding.step(callback_fn)

      // @ts-ignore
      ctx.embedding = embedding;
    } else if (e.data.messageType === 'abort') {
      // @ts-ignore
      const { embedding } = ctx;
      embedding.abort(callbackFnTerminated);
    }
    // else {
    //     // @ts-ignore
    //     const embedding = ctx.embedding;
    //     embedding.step(callback_fn)
    // }
  },
  false,
);

const helper = null;
export default helper;
