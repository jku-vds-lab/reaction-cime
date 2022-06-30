import { VisualizationSpec, VegaLite } from 'react-vega';
import { VegaLiteProps } from 'react-vega/lib/VegaLite';

const spec: VisualizationSpec = {
  "width": 50,
  "height": 50,
  "layer":[{
    "transform": [{
        "filter": "datum.selection == 'selection'"
      }],
    "mark": { "type": "area", "tooltip": true, "interpolate": "linear" },
    "encoding": {
      "x": { "field": "feature", "title": "feature", "type": "quantitative", "axis": null },
      "y": { "field": "density", "title": "density", "type": "quantitative", "axis": null },
      "color": { "field": "selection", "type": "nominal", "legend": null, "scale": { "range": ["(170,170,170,0)", "#007dad"] } }
    }
  },
  {
    "transform": [{
      "filter": "datum.selection == 'all'"
    }],
    "mark": { "type": "area", "tooltip": true, "interpolate": "linear" },
    "encoding": {
      "x": { "field": "feature", "title": "feature", "type": "quantitative", "axis": null },
      "y": { "field": "density", "title": "density", "type": "quantitative", "axis": null },
      "stroke": { "field": "selection", "type": "nominal", "legend": null, "scale": { "range": ["#000000", "#007dad"] } },
      "color": { "field": "selection", "type": "nominal", "legend": null }
    }
  }
  ],
  "data": { "name": 'values' },
};

export default function AreaChart (props: Omit<VegaLiteProps, 'spec'>) {
  return <VegaLite {...props} spec={spec} />;
}

// , "scale": {"domain": [0, 1]}
// rank != 1 and switching color scale range would be a workaround to make sure when there is only 1 datapoint it has the rank 1 colour
