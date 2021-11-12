import { Box, Typography } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import { setAggregateDatasetAction } from "../../State/AggregateDatasetDuck";
import { AppState } from "../../State/Store";
import { AggregatedDatasetDrop } from "./AggregatedDatasetDrop";

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  currentAggregation: state.currentAggregation,
  datasetEntries: state.datasetEntries
});

const mapDispatchToProps = (dispatch) => ({
  setDataset: dataset => dispatch(setAggregateDatasetAction(dataset)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  splitRef: any;
};

export const AggregationTabPanel = connector(
  ({setDataset, datasetEntries}: Props) => {
    console.log(datasetEntries)
    
    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <Typography variant="subtitle2" gutterBottom>
            Upload Aggregated Dataset
          </Typography>
        </Box>
        <AggregatedDatasetDrop
            onChange={(dataset) => {setDataset(dataset)}} />
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          asdf
        </Box>
      </div>
    );
  }
);

// function dispatch(arg0: { type: string; input: any }) {
//   throw new Error("Function not implemented.");
// }
