import * as React from 'react';
import { createRoot } from 'react-dom/client';
import { VisynAppProvider } from 'visyn_core/app';
import { ReactionCIMEApp } from './ReactionCIMEApp';
import './index.css';
import 'lineupjs/build/LineUpJS.css';

// Globally extend the clientConfig
declare module 'visyn_core/base' {
  export interface IClientConfig {
    publicVersion?: boolean;
  }
}

createRoot(document.getElementById('app')).render(
  <VisynAppProvider appName="CIME4R">
    <ReactionCIMEApp />
  </VisynAppProvider>,
);
