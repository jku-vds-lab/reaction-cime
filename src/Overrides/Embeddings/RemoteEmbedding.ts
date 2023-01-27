import type { IBaseProjection } from 'projection-space-explorer';
import { ReactionCIMEBackend } from '../../Backend/ReactionCIMEBackend';

export class RemoteEmbedding {
  private dataset: string = null;

  private params: object = null;

  private selected_feature_info: object = null;

  private embedding: IBaseProjection = null;

  private current_steps = 0;

  private abort_controller = null;

  private msg = '';

  private backend: ReactionCIMEBackend;

  constructor(backendUrl: string, dataset: string, params: object, init_embedding: IBaseProjection, selected_feature_info: object) {
    this.dataset = dataset;
    this.params = params;
    this.embedding = init_embedding;
    this.selected_feature_info = selected_feature_info;
    this.abort_controller = new AbortController();
    this.backend = new ReactionCIMEBackend(backendUrl, {});
  }

  // getUpdate(callback_fn:(emb, step, msg)=>void){
  //     callback_fn(this.embedding, this.current_steps, this.msg); //this.params["iterations"])
  // }

  jsonParseMultiple = (ss: string) => {
    ss = ss
      .split('\n')
      .map((l) => l.trim())
      .join('');
    let start = ss.indexOf('{');
    let open = 0;
    const res = [];
    for (let i = start; i < ss.length; i++) {
      if (ss[i] === '{' && (i < 2 || ss.slice(i - 2, i) !== '\\"')) {
        open++;
      } else if (ss[i] === '}' && (i < 2 || ss.slice(i - 2, i) !== '\\"')) {
        open--;
        if (open === 0) {
          res.push(JSON.parse(ss.substring(start, i + 1)));
          start = i + 1;
        }
      }
    }
    return res;
  };

  initializeFit(callback_fn: (emb, step, msg) => void) {
    // https://stackoverflow.com/questions/62121310/how-to-handle-streaming-data-using-fetch
    const readStream = async (reader) => {
      let done;
      let value;
      while (!done) {
        // eslint-disable-next-line no-await-in-loop
        ({ value, done } = await reader.read());
        if (done) {
          console.log('The stream is closed!');
        } else {
          try {
            const valueString = new TextDecoder().decode(value);
            // In case multiple JSON objects are received at once, we need to parse them all
            const splitValue = this.jsonParseMultiple(valueString);
            const resObj = splitValue[splitValue.length - 1];
            if (parseInt(resObj.step, 10)) this.current_steps = parseInt(resObj.step, 10);
            if (resObj.emb && resObj.emb.length > 0) this.embedding = resObj.emb;
            if (resObj.msg) this.msg = resObj.msg;
            callback_fn(this.embedding, this.current_steps, this.msg);
          } catch (error) {
            console.log(error);
            callback_fn(this.embedding, this.current_steps, this.msg);
          }
        }
      }
    };

    this.backend.project_dataset(this.dataset, this.params, this.selected_feature_info, this.abort_controller).then((response) => {
      const reader = response.body.getReader();
      readStream(reader);
    });
  }

  abort(callback) {
    // this.abort_controller.abort()
    this.backend
      .terminate_projection(this.dataset)
      .then((response) => {
        // if(response["msg"] === "ok"){
        this.abort_controller.abort();
        callback();
        // }
      })
      .catch(() => {
        this.abort_controller.abort();
        callback();
      });
  }
}
