import { defineConfig } from 'cypress';
import * as webpackConfig from 'visyn_scripts/config/webpack.config';

export default defineConfig({
  viewportHeight: 1080,
  viewportWidth: 1920,
  defaultCommandTimeout: 10000,
  e2e: {
    baseUrl: 'http://localhost:8080',
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    setupNodeEvents(on, config) {},
  },
  component: {
    devServer: {
      framework: 'react',
      bundler: 'webpack',
      webpackConfig: webpackConfig({ workspace_mode: 'single' }, { mode: 'production' }),
    },
  },
});
