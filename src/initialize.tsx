import * as React from 'react';
import ReactDOM from 'react-dom';
import { Anchor } from '@mantine/core';
import { VisynApp, VisynAppProvider, VisynHeader } from 'visyn_core/app';
import { ReactionCIMEApp } from './ReactionCIMEApp';
import './index.css';
import 'lineupjs/build/LineUpJS.css';
import vdsLogo from './assets/jku-vds-lab-logo.svg';
import bayerLogo from './assets/bayer_logo.svg';

ReactDOM.render(
  <VisynAppProvider appName="CIME4R">
    <VisynApp
      header={
        <VisynHeader
          components={{
            beforeRight: (
              <>
                <Anchor
                  href="https://www.bayer.com/"
                  rel="noreferrer"
                  target="_blank"
                  sx={{
                    // Center the image
                    display: 'flex',
                    alignItems: 'center',
                  }}
                >
                  <img src={bayerLogo} alt="Bayer logo" style={{ height: '32px', transform: 'scale(1.1)' }} />
                </Anchor>
                <Anchor
                  href="https://jku-vds-lab.at/"
                  rel="noreferrer"
                  target="_blank"
                  sx={{
                    // Center the image
                    display: 'flex',
                    alignItems: 'center',
                  }}
                >
                  <img src={vdsLogo} alt="JKU VDS Lab logo" style={{ height: '24px' }} />
                </Anchor>
              </>
            ),
          }}
        />
      }
    >
      <ReactionCIMEApp />
    </VisynApp>
  </VisynAppProvider>,
  document.getElementById('app'),
);
