import * as React from 'react';
import ReactDOM from 'react-dom';
import { VisynAppProvider } from 'visyn_core/app';
import { ReactionCIMEApp } from './ReactionCIMEApp';
import './index.css';
import 'lineupjs/build/LineUpJS.css';

ReactDOM.render(
  <VisynAppProvider
    appName="CIME4R"
    defaultClientConfig={{
      // The visyn app will load a clientConfig.json from the root of the webservice, and override those properties here. See the Dockerfile for more information.
      publicVersion: false,
    }}
  >
    <ReactionCIMEApp />
  </VisynAppProvider>,
  document.getElementById('app'),
);
