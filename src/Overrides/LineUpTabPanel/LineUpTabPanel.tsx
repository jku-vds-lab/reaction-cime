import { Typography, Box, Button, FormControlLabel, Switch } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import GetAppIcon from '@mui/icons-material/GetApp';
import React from 'react';
import { DetailViewActions } from 'projection-space-explorer';
import { setLineUpInputFilter } from '../../State/LineUpInputDuck';
import { AppState } from '../../State/Store';
import { downloadImpl } from '../../Utility/Utils';

const mapStateToProps = (state: AppState) => ({
  currentAggregation: state.currentAggregation,
  lineUpInput: state.lineUpInput,
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

export const LineUpTabPanel = connector(({ setDetailVisibility, setLineUpInput_filter, lineUpInput, currentAggregation, splitRef }: Props) => {
  // const handleChange = (_, value) => {};

  const onLoad = (filter) => {
    setDetailVisibility(true);
    setLineUpInput_filter(filter);

    const currSizes = splitRef.current.split.getSizes();
    if (currSizes[1] < 2) {
      splitRef.current.split.setSizes([currSizes[0], 70]);
    }
  };

  // https://stackoverflow.com/questions/31214677/download-a-reactjs-object-as-a-file
  // const downloadImpl = (data: string, name: string, mimetype: string) => {
  //   var b = new Blob([data], { type: mimetype });
  //   var csvURL = window.URL.createObjectURL(b);
  //   let tempLink = document.createElement("a");
  //   tempLink.href = csvURL;
  //   tempLink.setAttribute("download", name);
  //   tempLink.click();
  // };

  const exportCSV = () => {
    if (lineUpInput.lineup && lineUpInput.lineup.data) {
      // exports all data that is currently shown in the table -> filters and sorts are applied! also annotations are included
      lineUpInput
        .lineup!.data.exportTable(lineUpInput.lineup!.data.getRankings()[0], {
          separator: ',',
        })
        .then((x) => downloadImpl(x, `lineup-export.csv`, 'application/csv'));
    }
  };

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
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Typography variant="subtitle2" gutterBottom>
          LineUp Settings
        </Typography>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        {/* <FormControlLabel
                control={<Switch checked={true} onChange={handleChange} name="checkedA" />}
                label="External Selection Summary"
            /> */}

        <Button fullWidth style={{ marginRight: 2 }} variant="outlined" onClick={() => onLoad({ reset: true })}>
          Load All
        </Button>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        <Button fullWidth variant="outlined" onClick={() => onLoad({ selection: currentAggregation.aggregation })}>
          Load Selection
        </Button>
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
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
          label="Show Cell Values"
        />
      </Box>
      <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
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
      </Box>
    </div>
  );
});
