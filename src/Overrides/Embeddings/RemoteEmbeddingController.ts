import { Dataset, EmbeddingController, IBaseProjection } from "projection-space-explorer";
import remoteWorker from "./remote.worker";

export class RemoteEmbeddingController extends EmbeddingController {
    targetBounds: any
    private embedding_method:string;

    constructor(embedding_method:string){
        super()
        this.embedding_method = embedding_method
    }
    
    init(dataset: Dataset, selection: any, params: any, workspace: IBaseProjection) {
        const selected_feature_info = {}
        selection.forEach(feature => {
            if(feature.checked){
                let feature_info = {};
                feature_info["normalize"] = feature.normalized;
                feature_info["featureType"] = dataset.columns[feature.name].featureType.toString();
                feature_info["range"] = dataset.columns[feature.name].range;
                feature_info["distinct"] = dataset.columns[feature.name].distinct;
                selected_feature_info[feature.name] = feature_info;
            }
        });
        params["embedding_method"] = this.embedding_method;

        this.worker = new remoteWorker()
        this.worker.postMessage({
            messageType: 'init',
            init_coordinates: workspace.map((v, i) => {
                return {x: v.x, y: v.y}
            }),
            path: dataset.info.path,
            params: params,
            selected_feature_info: selected_feature_info
        })

        this.worker.addEventListener('message', (e) => {
            var Y = e.data
            this.stepper(Y)
            this.notifier()
        }, false);
    }

    step() {
        this.worker.postMessage({
            messageType: 'step'
        })
    }
}