import { Box, Typography } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { AggregatedDatasetDrop } from "./AggregatedDatasetDrop";

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  currentAggregation: state.currentAggregation,
});

const mapDispatchToProps = (dispatch) => ({
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  splitRef: any;
};

export const AggregationTabPanel = connector(
  ({}: Props) => {
    
    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <Typography variant="subtitle2" gutterBottom>
            Upload Aggregated Dataset
          </Typography>
        </Box>
        <AggregatedDatasetDrop
            onChange={() => {console.log("test")}} />
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          ASDF
        </Box>
      </div>
    );
  }
);

// function dispatch(arg0: { type: string; input: any }) {
//   throw new Error("Function not implemented.");
// }
