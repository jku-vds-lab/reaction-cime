import { EmbeddingController } from 'projection-space-explorer';

export class DummyEmbeddingController extends EmbeddingController {
  init() {
    console.log('init');
  }

  override terminate() {
    console.log('terminate');
  }

  override supportsPause() {
    console.log('supportsPause');
    return false;
  }
}
