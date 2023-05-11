import type { IBaseProjection } from 'projection-space-explorer';
import { ReactionCIMEBackend } from '../../Backend/ReactionCIMEBackend';

export class RemoteEmbedding {
  private current_steps = 0;

  private abort_controller = null;

  private msg = '';

  private backend: ReactionCIMEBackend;

  constructor(
    backendUrl: string,
    protected dataset: string,
    protected params: object,
    protected embedding: IBaseProjection,
    protected selected_feature_info: object,
    protected ids: number[],
  ) {
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
    if (start === -1) {
      return { res, rest: ss };
    }
    for (let i = start; i < ss.length; i++) {
      if (ss[i] === '{' && (i < 2 || ss.slice(i - 2, i) !== '\\"')) {
        open++;
      } else if (ss[i] === '}' && (i < 2 || ss.slice(i - 2, i) !== '\\"')) {
        open--;
        if (open === 0) {
          try {
            res.push(JSON.parse(ss.substring(start, i + 1)));
            start = i + 1;
          } catch (e) {
            console.error(`Could not properly parse JSON object from '${ss}', ignoring it`, e);
            // ignore
          }
        }
      }
    }
    return { parsed: res, rest: ss.substring(start).trim() };
  };

  initializeFit(callback_fn: (emb, step, msg) => void) {
    this.backend
      .project_dataset(this.dataset, this.params, this.selected_feature_info, this.ids, this.abort_controller)
      .then(async (response) => {
        // https://stackoverflow.com/questions/62121310/how-to-handle-streaming-data-using-fetch
        const reader = response.body.getReader();
        let done;
        let value;
        let previousValue = '';
        /**
         * True if the embedding has been set at least once, used to show the user an error if they should reload the page.
         */
        let hasSetEmbedding = false;
        while (!done) {
          // eslint-disable-next-line no-await-in-loop
          ({ value, done } = await reader.read());
          if (done) {
            if (!hasSetEmbedding) {
              try {
                // The embedding was not set automatically, so we fetch it here
                // eslint-disable-next-line no-await-in-loop
                const coordinates = await this.backend.get_coordinates(this.dataset, this.ids);
                callback_fn(coordinates, this.current_steps, this.msg);
              } catch (e) {
                console.error('Error while fetching embedding', e);
                callback_fn(this.embedding, this.current_steps, 'Embedding could not be loaded, please reload the application to see the new embedding.');
              }
            }
          } else {
            try {
              const valueString = new TextDecoder().decode(value);
              // In case multiple JSON objects are received at once, we need to parse them all
              const { parsed: splitValue, rest } = this.jsonParseMultiple(previousValue + valueString);
              if (splitValue.length === 0) continue;
              previousValue = rest;
              const resObj = splitValue[splitValue.length - 1];
              if (parseInt(resObj.step, 10)) this.current_steps = parseInt(resObj.step, 10);
              if (resObj.emb && resObj.emb.length > 0) {
                hasSetEmbedding = true;
                this.embedding = resObj.emb;
              }
              if (resObj.msg) this.msg = resObj.msg;
              callback_fn(this.embedding, this.current_steps, this.msg);
            } catch (error) {
              console.error('Error while projecting', error);
              callback_fn(this.embedding, this.current_steps, this.msg);
            }
          }
        }
      })
      .catch((e) => {
        console.error('Error while projecting', e);
        callback_fn(this.embedding, this.current_steps, `Error while projecting: ${e}`);
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
