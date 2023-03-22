import { Typography, Box, Button, FormControlLabel, Switch, Tooltip } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import React from 'react';
import { DetailViewActions } from 'projection-space-explorer';
import { setLineUpInputFilter } from '../../State/LineUpInputDuck';
import { AppState } from '../../State/Store';

const mapStateToProps = (state: AppState) => ({
  currentAggregation: state.currentAggregation,
  lineUpInput: state.lineUpInput,
  globalLabels: state.globalLabels,
});

const mapDispatchToProps = (dispatch) => ({
  setDetailVisibility: (value) => dispatch(DetailViewActions.setDetailVisibility(value)),
  setLineUpInput_filter: (value) => dispatch(setLineUpInputFilter(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  splitRef: any;
};

export const LineUpTabPanel = connector(({ setDetailVisibility, setLineUpInput_filter, lineUpInput, currentAggregation, splitRef, globalLabels }: Props) => {
  // const handleChange = (_, value) => {};

  const onLoad = (filter) => {
    setDetailVisibility(true);
    setLineUpInput_filter(filter);

    const currSizes = splitRef.current.split.getSizes();
    if (currSizes[1] < 2) {
      splitRef.current.split.setSizes([currSizes[0], 70]);
    }
  };

  // const exportCSV = () => {
  //   if (lineUpInput.lineup && lineUpInput.lineup.data) {
  //     // exports all data that is currently shown in the table -> filters and sorts are applied! also annotations are included
  //     lineUpInput
  //       .lineup!.data.exportTable(lineUpInput.lineup!.data.getRankings()[0], {
  //         separator: ',',
  //       })
  //       .then((x) => downloadImpl(x, `lineup-export.csv`, 'application/csv'));
  //   }
  // };

  const [cellValueVis, setCellValueVis] = React.useState(false);

  React.useEffect(() => {
    let style = document.getElementById('cell_value_vis');
    if (!style) {
      style = document.createElement('style');
      style.setAttribute('id', 'cell_value_vis');
      style.setAttribute('type', 'text/css');
      const head = document.head || document.getElementsByTagName('head')[0];
      head.appendChild(style);
    }

    const css = cellValueVis ? '.lu-hover-only { visibility: visible; }' : '.lu-hover-only { visibility: hidden; }';
    // @ts-ignore
    if (style.styleSheet) {
      // This is required for IE8 and below.
      // @ts-ignore
      style.styleSheet.cssText = css;
    } else {
      style.innerHTML = '';
      style.appendChild(document.createTextNode(css));
    }
  }, [cellValueVis]);

  const toggleVis = () => {
    setCellValueVis(() => {
      return !cellValueVis;
    });
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      <Box paddingX={2} paddingTop={2} paddingBottom={1}>
        <Typography variant="subtitle2" gutterBottom>
          LineUp settings
        </Typography>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Tooltip placement="right" title={<Typography variant="subtitle2">Show all {globalLabels.itemLabelPlural} in the LineUp view</Typography>}>
          <Button fullWidth style={{ marginRight: 2 }} variant="outlined" onClick={() => onLoad({ reset: true })}>
            View all {globalLabels.itemLabelPlural}
          </Button>
        </Tooltip>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Tooltip
          placement="right"
          title={<Typography variant="subtitle2">Only show the current selection of {globalLabels.itemLabelPlural} in the LineUp view</Typography>}
        >
          <span>
            <Button
              disabled={currentAggregation?.aggregation.length <= 0}
              fullWidth
              variant="outlined"
              onClick={() => onLoad({ selection: currentAggregation.aggregation })}
            >
              View selected {globalLabels.itemLabelPlural}
            </Button>
          </span>
        </Tooltip>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Tooltip placement="right" title={<Typography variant="subtitle2">If activated, LinUp shows values in the cells for numeric features.</Typography>}>
          <FormControlLabel
            control={
              <Switch
                color="primary"
                value={cellValueVis}
                onChange={(event) => {
                  toggleVis();
                }}
              />
            }
            label="Show cell values"
          />
        </Tooltip>
      </Box>
      {/* TODO: if we need this, debug before activating again */}
      {/* <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Button
          fullWidth
          variant="outlined"
          onClick={() => {
            exportCSV();
          }}
        >
          <GetAppIcon />
          &nbsp;Export CSV
        </Button>
      </Box> */}
    </div>
  );
});
