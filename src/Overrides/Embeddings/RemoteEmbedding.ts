import { IBaseProjection } from "projection-space-explorer";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";


export class RemoteEmbedding {
    private dataset:string = null;
    private params:object = null;
    private selected_feature_info:object = null;
    private embedding:IBaseProjection = null;
    private current_steps = 0;
    private abort_controller = null;
    private msg = "";

    constructor(dataset:string, params:object, init_embedding:IBaseProjection, selected_feature_info:object) {
        this.dataset = dataset;
        this.params = params;
        this.embedding = init_embedding;
        this.selected_feature_info = selected_feature_info;
        this.abort_controller = new AbortController();
        
    }

    // getUpdate(callback_fn:(emb, step, msg)=>void){
    //     callback_fn(this.embedding, this.current_steps, this.msg); //this.params["iterations"])
    // }

    initializeFit(callback_fn:(emb, step, msg)=>void){
        // https://stackoverflow.com/questions/62121310/how-to-handle-streaming-data-using-fetch
        const read_stream = async(reader) => {
            let done, value;
            while (!done) {
              ({ value, done } = await reader.read());
              if (done) {
                console.log("The stream is closed!");
              } else {
                  try{
                    const res_obj = JSON.parse(new TextDecoder().decode(value))
                    if(parseInt(res_obj["step"]))
                        this.current_steps = parseInt(res_obj["step"]);
                    if(res_obj["emb"] && res_obj["emb"].length > 0)
                        this.embedding = res_obj["emb"];
                    if(res_obj["msg"])
                        this.msg = res_obj["msg"]
                    callback_fn(this.embedding, this.current_steps, this.msg);
                  } catch (error){
                      console.log(error);
                      callback_fn(this.embedding, this.current_steps, this.msg);
                  }
              }
            }
          };

        ReactionCIMEBackendFromEnv.project_dataset(this.dataset, this.params, this.selected_feature_info, this.abort_controller).then((response)=>{
            const reader = response.body.getReader()
            read_stream(reader)
        })
    }

    abort(callback){
        // this.abort_controller.abort()
        ReactionCIMEBackendFromEnv.terminate_projection(this.dataset).then((response) => {
            // if(response["msg"] === "ok"){
                this.abort_controller.abort()
                callback()
            // }
        }).catch(() => {
            this.abort_controller.abort()
            callback()
        })
    }

    // step(callback_fn:(emb)=>void){
    //     // TODO: implement step logic
    //     this.getUpdate(callback_fn)
    // }
}