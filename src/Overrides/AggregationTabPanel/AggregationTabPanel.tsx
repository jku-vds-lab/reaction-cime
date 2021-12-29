import { Box } from "@mui/material";
import { CategoryOptionsAPI, SelectFeatureComponent } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import { useEffect } from "react";

const mapStateToProps = (state: AppState) => ({
  aggregateColor: state.aggregateColor,
  poiDataset: state.dataset,
  workspace: state.projections.workspace
});

const mapDispatchToProps = (dispatch) => ({
  setAggregateColor: aggregateColor => dispatch(setAggregateColor(aggregateColor)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};

export const AggregationTabPanel = connector(({setAggregateColor, aggregateColor, poiDataset, workspace}: Props) => {
    const categoryOptions = poiDataset?.categories;

    useEffect(() => {
      // reset aggregate color to hide the aggregated dataset in the background
      setAggregateColor({key: "None", name: "None"})
    // eslint-disable-next-line
    }, [workspace]) // this is triggered during the embedding

    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
        { //TODO: check, if it makes sense to also include categorical values, or if it is ok to only use numerical values (like for "size")
          categoryOptions != null && CategoryOptionsAPI.hasCategory(categoryOptions, "size") ?
              <SelectFeatureComponent label={"color"} default_val={aggregateColor} categoryOptions={CategoryOptionsAPI.getCategory(categoryOptions, "size")} onChange={(newValue) => {
                  var attribute = null
                  if (newValue && newValue !== "") {
                    attribute = CategoryOptionsAPI.getCategory(categoryOptions, "size").attributes.filter(a => a.key === newValue)[0]
                  }
                  if (attribute === null || attribute === undefined){
                    attribute = {key: "None", name: "None"}
                  }
                  setAggregateColor(attribute);

              }}></SelectFeatureComponent>
              :
              <div></div>
        }
        </Box>
      </div>
    );
  }
);

