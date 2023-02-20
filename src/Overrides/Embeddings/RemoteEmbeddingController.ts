import { Dataset, EmbeddingController, IProjection } from 'projection-space-explorer';
// @ts-ignore
import RemoteWorker from 'worker-loader?inline=no-fallback!./remote.worker';
import { ReactionCIMEBackendFromEnv } from '../../Backend';

export class RemoteEmbeddingController extends EmbeddingController {
  targetBounds: any;

  private embedding_method: string;

  private callback_fn?: (msg: string) => void;

  constructor(embedding_method: string, callback_fn?: (msg: string) => void) {
    // embedding_methods currently implemented in the backend: "umap", "tsne", "pca"; it defaults to "pca"
    super();
    this.embedding_method = embedding_method;
    this.callback_fn = callback_fn;
  }

  init(dataset: Dataset, selection: any, params: any, workspace: IProjection) {
    if (this.callback_fn) this.callback_fn('init');
    const selectedFeatureInfo = {};
    selection.forEach((feature) => {
      if (feature.checked) {
        const featureInfo = {
          normalize: feature.normalized,
          featureType: dataset.columns[feature.name].featureType.toString(),
          range: dataset.columns[feature.name].range,
          distinct: dataset.columns[feature.name].distinct,
          useWeight: feature.useWeight,
          weight: feature.weight,
        };
        selectedFeatureInfo[feature.name] = featureInfo;
      }
    });
    params.embedding_method = this.embedding_method;
    this.worker = new RemoteWorker();
    this.worker.postMessage({
      messageType: 'init',
      backendUrl: ReactionCIMEBackendFromEnv.baseUrl,
      init_coordinates: workspace.positions.map((v, i) => {
        return { x: v.x, y: v.y };
      }),
      path: dataset.info.path,
      params,
      selected_feature_info: selectedFeatureInfo,
    });

    this.worker.addEventListener(
      'message',
      (e) => {
        if (e.data.messageType === 'step') {
          this.stepper(e.data.embedding);
          // @ts-ignore
          this.notifier(e.data.step, e.data.msg);
        } else if (e.data.messageType === 'terminated') {
          this.worker.terminate();
        }
      },
      false,
    );
  }

  terminate() {
    this.worker.postMessage({
      messageType: 'abort',
    });
  }

  step() {
    // step is not controlled by the ProjectionControlCard anymore
    // this.worker.postMessage({
    //     messageType: 'step'
    // })
  }
}
