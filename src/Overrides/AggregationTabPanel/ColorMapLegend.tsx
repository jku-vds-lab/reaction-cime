import * as React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { Button, InputLabel, MenuItem, Select, Tooltip, Typography } from '@mui/material';
import * as vsup from 'vsup';
import * as d3 from 'd3v5';
import { InfoOutlined } from '@mui/icons-material';
import { Box } from '@mui/system';
import { D3_CONTINUOUS_COLOR_SCALE_LIST, AggregateActions } from '../../State/AggregateSettingsDuck';
import { AppState } from '../../State/Store';
import { PSE_BLUE } from '../../Utility/Utils';

const mapStateToProps = (state: AppState) => ({
  aggregateSettings: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings,
  selectAttribute: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttribute,
  activeId: state.multiples.active,
});

const mapDispatchToProps = (dispatch) => ({
  setAggregateColorScale: (value) => dispatch(AggregateActions.setAggregateColorScale(value)),
  toggleUseVSUP: () => dispatch(AggregateActions.toggleUseVSUP()),
  addValueFilter: (value) => dispatch(AggregateActions.addValueFilter(value)),
  removeValueFilter: (value) => dispatch(AggregateActions.removeValueFilter(value)),
  clearValueFilter: () => dispatch(AggregateActions.clearValueFilter()),
  setAggregateColorMapScale: (scale) => dispatch(AggregateActions.setAggregateColorMapScale(scale)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

export const ColorMapLegend = connector(
  ({
    setAggregateColorScale,
    selectAttribute,
    toggleUseVSUP,
    addValueFilter,
    removeValueFilter,
    clearValueFilter,
    aggregateSettings,
    setAggregateColorMapScale,
    activeId,
  }: Props) => {
    if (selectAttribute == null || selectAttribute === 'None') {
      return null;
    }

    const discreteSteps = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

    const handleChange = (event) => {
      setAggregateColorScale(event.target.value);
    };

    // add colormap legend
    const gRef = React.useRef();
    const svgRef = React.useRef();
    const [colorSections, setColorSections] = React.useState([]);
    const [legend, setLegend] = React.useState(null);

    const clearFilter = (color_sections) => {
      clearValueFilter();

      color_sections.attr('stroke', null);
      color_sections.attr('stroke-width', '0px');
    };

    React.useEffect(() => {
      if (aggregateSettings.advancedSettings.valueRange != null) {
        // visit https://github.com/uwdata/vsup for more info about VSUP vs bivariate colorscale encoding
        let vDom = [aggregateSettings.advancedSettings.valueRange.min, aggregateSettings.advancedSettings.valueRange.max];
        let colLegend;
        let scale;
        let quantization;

        // simple encoding
        if (aggregateSettings.advancedSettings.uncertaintyRange == null) {
          scale = d3
            .scaleQuantize()
            // @ts-ignore
            .domain(vDom)
            .range(d3.quantize(d3[aggregateSettings.colormapSettings.colorscale], 8));
          vDom = scale.domain();
          colLegend = vsup.legend.simpleLegend(scale).title(aggregateSettings.colormapSettings.aggregateColor.value_col);
        } else {
          // bivariate encoding
          const uDom = [aggregateSettings.advancedSettings.uncertaintyRange.min, aggregateSettings.advancedSettings.uncertaintyRange.max];
          if (aggregateSettings.colormapSettings.useVSUP) {
            quantization = vsup.quantization().branching(2).layers(4).valueDomain(vDom).uncertaintyDomain(uDom);
            scale = vsup.scale().quantize(quantization).range(d3[aggregateSettings.colormapSettings.colorscale]);
            colLegend = vsup.legend.arcmapLegend(scale);
          } else {
            quantization = vsup.squareQuantization(4).valueDomain(vDom).uncertaintyDomain(uDom);
            scale = vsup.scale().quantize(quantization).range(d3[aggregateSettings.colormapSettings.colorscale]);
            colLegend = vsup.legend.heatmapLegend(scale);
          }
          colLegend
            .vtitle(aggregateSettings.colormapSettings.aggregateColor.value_col)
            .utitle(aggregateSettings.colormapSettings.aggregateColor.uncertainty_col);
        }

        // TODO: for extremely small values, we could add some scaling function that scales the values by some e-5 and add "e-5" to the label --> the ticks would then for example be "0.01" and the title "columnname e-5"
        // automatically obtain a suitable precision value
        const p = d3.precisionFixed((vDom[1] - vDom[0]) / 8);
        if (p > 2) {
          colLegend.format('.2');
        } else if (p > 0) {
          colLegend.format(`.${p}r`);
        }

        setLegend(() => colLegend);
        setAggregateColorMapScale(scale);
      }
      // eslint-disable-next-line
    }, [
      activeId,
      aggregateSettings.colormapSettings.useVSUP,
      aggregateSettings.advancedSettings.valueRange,
      aggregateSettings.advancedSettings.uncertaintyRange,
      aggregateSettings.colormapSettings.colorscale,
      aggregateSettings.colormapSettings.aggregateColor,
    ]);

    React.useEffect(() => {
      if (legend != null && gRef.current && svgRef.current) {
        const relWidth = 250;
        legend.size(relWidth);

        let paddingLeft;
        let paddingTop;
        let paddingRight;
        let paddingBottom;
        let relHeight;
        if (legend.height == null) {
          relHeight = relWidth;
          paddingLeft = relWidth * 0.1;
          paddingTop = paddingLeft * 2;
          paddingRight = paddingLeft * 3;
          paddingBottom = paddingLeft * 3;
        } else {
          relHeight = legend.height();
          paddingLeft = relWidth * 0.05;
          paddingTop = 0;
          paddingRight = paddingLeft * 2;
          paddingBottom = relHeight * 1.2;
        }

        const gElement = d3.select(gRef.current);
        gElement.html(''); // clear g element
        gElement.call(legend); // draw color legend

        const svgElement = d3.select(svgRef.current);
        svgElement.attr('viewBox', `-${paddingLeft} -${paddingTop} ${relWidth + paddingRight} ${relHeight + paddingBottom}`);

        // --- add interaction with legend
        const legendContainer = svgElement.select('.legend > g:last-child');
        let newColorSections = legendContainer.selectAll('path');
        if (newColorSections.nodes().length <= 0) {
          // if there are no path elements, we look for rect elements
          newColorSections = gElement.selectAll('rect');
        }
        setColorSections(newColorSections.nodes());

        newColorSections.on('mouseover', (d, i) => {
          // TODO: add hover interaction --> highlight areas in background somehow (e.g. have a component that draws a border around the selected areas)
          legendContainer.select('.hover_clone').remove();
          const hoverEl = d3.select(newColorSections.nodes()[i]).clone();
          hoverEl.attr('class', 'hover_clone');
          hoverEl.attr('fill', PSE_BLUE);

          // remove temporary hover element when we leave it
          hoverEl.on('mouseout', () => {
            hoverEl.remove();
          });

          // when clicking, we want to select the underlying section
          hoverEl.on('click', () => {
            // color_sections.attr("stroke", "none")
            const curNode = d3.select(newColorSections.nodes()[i]);
            if (curNode.attr('stroke') == null) {
              curNode.attr('stroke', PSE_BLUE);
              curNode.attr('stroke-width', '3px');
              addValueFilter(curNode.attr('fill'));
            } else {
              curNode.attr('stroke', null);
              curNode.attr('stroke-width', '0px');
              removeValueFilter(curNode.attr('fill'));
            }
          });
        });

        // svgElement.on("mouseout", () => {
        // legend_container.select(".hover_clone").remove();
        // });
      }
      // eslint-disable-next-line
    }, [legend, gRef, svgRef]);

    return (
      <>
        <Typography variant="body2" color="textSecondary">
          Choose colormap
        </Typography>
        <Select labelId="colorscale-select-label" id="colorscale-select" value={aggregateSettings.colormapSettings.colorscale} onChange={handleChange}>
          {D3_CONTINUOUS_COLOR_SCALE_LIST.map((colorscale) => (
            <MenuItem key={colorscale} value={colorscale} title={colorscale}>
              <div
                style={{
                  width: '100%',
                  minWidth: '13rem',
                  height: '1rem',
                  backgroundImage: `linear-gradient(to right, ${discreteSteps.map((step) => d3[colorscale](step)).join(',')})`,
                }}
              />
            </MenuItem>
          ))}
        </Select>
        <svg ref={svgRef} style={{ cursor: 'pointer' }}>
          <g ref={gRef} />
        </svg>
        {legend?.height == null && ( // only the "simple" legend has a hight attribute
          <Box>
            <Typography variant="body2" color="textSecondary" gutterBottom>
              The chosen feature consists of two variables. You can choose bi-variate color mapping or VSUP.{' '}
              <Tooltip
                placement="right"
                title={
                  <Typography variant="subtitle2">
                    Both color mapping methods can be used for encoding two variables. In bi-variate color mapping, the two variables are linearly encoded with
                    hue and saturation. For value-suppressing uncertainty palettes (VSUP), one of the values is considered a measure of &quot;uncertainty&quot;.
                    With increasing uncertainty the separation of values becomes less important.
                  </Typography>
                }
              >
                <InfoOutlined fontSize="inherit" style={{ color: 'grey' }} />
              </Tooltip>
            </Typography>
            <Tooltip
              placement="right"
              title={<Typography variant="subtitle2">Switch to {aggregateSettings.colormapSettings.useVSUP ? 'bi-variate mapping' : 'VSUP'}</Typography>}
            >
              <Button
                variant="outlined"
                onClick={() => {
                  toggleUseVSUP();
                }}
              >
                Switch encoding
              </Button>
            </Tooltip>
          </Box>
        )}
        {aggregateSettings.colormapSettings.valueFilter != null && aggregateSettings.colormapSettings.valueFilter.length > 0 && (
          <Button
            variant="outlined"
            onClick={() => {
              clearFilter(colorSections);
            }}
          >
            Clear filter
          </Button>
        )}
      </>
    );
  },
);
