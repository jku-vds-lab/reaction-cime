/* global require, module, __dirname */
// const { override, fixBabelImports, addWebpackAlias } = require('customize-cra')
const path = require('path')

// module.exports = override(
//   addWebpackAlias({
//     'react': path.resolve('./node_modules/react'),
//     'react-dom': path.resolve('./node_modules/react-dom'),
//     'react-redux': path.resolve('./node_modules/react-redux'),
//     'redux': path.resolve('./node_modules/redux'),
//     '@mui/material': path.resolve('./node_modules/@mui/material'),
//     '@emotion/styled': path.resolve('./node_modules/@emotion/styled'),
//     '@emotion/react': path.resolve('./node_modules/@emotion/react'),
//     '@mui/styles': path.resolve('./node_modules/@mui/styles'),
//   })
// )


module.exports = function override(config, env) {
  config.resolve.alias = {
        'react': path.resolve('./node_modules/react'),
        'react-dom': path.resolve('./node_modules/react-dom'),
        'react-redux': path.resolve('./node_modules/react-redux'),
        'redux': path.resolve('./node_modules/redux'),
        '@mui/material': path.resolve('./node_modules/@mui/material'),
        '@emotion/styled': path.resolve('./node_modules/@emotion/styled'),
        '@emotion/react': path.resolve('./node_modules/@emotion/react'),
        '@mui/styles': path.resolve('./node_modules/@mui/styles'),
      }

      // need to add it to beginning of array i.e. *before* the `ts-loader`
  // config.module.rules.unshift({
    config.module.rules.push({
      test: /\.worker\.ts$/,
      loader: 'worker-loader',
      options: { inline: "no-fallback" }
    })
    config.module.rules.push({
        test: /\.svg$/,
        use: ['@svgr/webpack'],
    })

  return config;
}

