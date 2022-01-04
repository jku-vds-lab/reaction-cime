import { Dataset, EmbeddingController, IBaseProjection } from "projection-space-explorer";
import remoteWorker from "./remote.worker";

export class RemoteEmbeddingController extends EmbeddingController {
    targetBounds: any;
    private embedding_method:String;

    constructor(embedding_method:String){ 
        // embedding_methods currently implemented in the backend: "umap", "tsne", "pca"; it defaults to "pca"
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
            if(e.data.messageType === "step"){
                this.stepper(e.data["embedding"])
                this.notifier(e.data["step"], e.data["msg"])
            }else if(e.data.messageType === "terminated"){
                this.worker.terminate();
            }
        }, false);
    }

    terminate(){
        this.worker.postMessage({
            messageType: 'abort'
        });
        
    }

    step() {
        // step is not controlled by the ProjectionControlCard anymore
        // this.worker.postMessage({
        //     messageType: 'step'
        // })
    }
}