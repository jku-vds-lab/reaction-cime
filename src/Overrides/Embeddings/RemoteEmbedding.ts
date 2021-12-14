import { IBaseProjection } from "projection-space-explorer";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";


export class RemoteEmbedding {
    private dataset:string = null;
    private params:object = null;
    private embedding:IBaseProjection = null;
    private selected_feature_info:object = null;

    constructor(dataset:string, params:object, init_embedding:IBaseProjection, selected_feature_info:object) {
        this.dataset = dataset;
        this.params = params;
        this.embedding = init_embedding;
        this.selected_feature_info = selected_feature_info;
    }

    getEmbedding(callback_fn:(emb)=>void){
        callback_fn(this.embedding)
    }

    initializeFit(callback_fn:(emb)=>void){
        ReactionCIMEBackendFromEnv.project_dataset(this.dataset, this.params, this.selected_feature_info).then((response)=>{
            if(response.done === "true"){
                this.embedding = response.embedding;
            }
            this.getEmbedding(callback_fn);
        })
    }

    step(callback_fn:(emb)=>void){
        // TODO: implement step logic
        this.getEmbedding(callback_fn)
    }
}