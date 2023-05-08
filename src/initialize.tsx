import * as React from 'react';
import ReactDOM from 'react-dom';
import { VisynAppProvider } from 'visyn_core/app';
import { ReactionCIMEApp } from './ReactionCIMEApp';
import './index.css';
import 'lineupjs/build/LineUpJS.css';

// Globally extend the clientConfig
declare module 'visyn_core' {
  export interface IClientConfig {
    publicVersion?: boolean;
  }
}

ReactDOM.render(
  <VisynAppProvider appName="CIME4R">
    <ReactionCIMEApp />
  </VisynAppProvider>,
  document.getElementById('app'),
);
