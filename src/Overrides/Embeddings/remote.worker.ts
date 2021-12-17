// import "regenerator-runtime/runtime";

import { RemoteEmbedding } from "./RemoteEmbedding";


/* eslint-disable-next-line no-restricted-globals */
const ctx = self;


const callback_fn_step = (emb, step, msg) => {ctx.postMessage({"messageType":"step", "embedding": emb, "step": step, "msg": msg})};
const callback_fn_terminated = () => {ctx.postMessage({"messageType":"terminated"})};

ctx.addEventListener('message', function (e) {
    if (e.data.messageType === 'init') {
        const embedding = new RemoteEmbedding(e.data.path, e.data.params, e.data.init_coordinates, e.data.selected_feature_info);
        embedding.initializeFit(callback_fn_step)
        // embedding.step(callback_fn)

        // @ts-ignore
        ctx.embedding = embedding;
    } else if (e.data.messageType === "abort"){
        // @ts-ignore
        const embedding = ctx.embedding;
        embedding.abort(callback_fn_terminated)
    }
    // else {
    //     // @ts-ignore
    //     const embedding = ctx.embedding;
    //     embedding.step(callback_fn)
    // }
}, false);


export default null;