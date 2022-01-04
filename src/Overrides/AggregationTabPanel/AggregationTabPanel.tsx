import { TextField, Box } from "@mui/material";
import { CategoryOptionsAPI, SelectFeatureComponent } from "projection-space-explorer";
// import { CategoryOptionsAPI, SelectFeatureComponent, useCancellablePromise } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import { setDatasetAction } from "projection-space-explorer/dist/components/Ducks/DatasetDuck"
import React from "react";


const mapStateToProps = (state: AppState) => ({
  aggregateColor: state.aggregateColor,
  poiDataset: state.dataset,
  workspace: state.projections.workspace
});

const mapDispatchToProps = (dispatch) => ({
  setDataset: dataset => dispatch(setDatasetAction(dataset)),
  setAggregateColor: aggregateColor => dispatch(setAggregateColor(aggregateColor)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};


export const AggregationTabPanel = connector(
  ({
    setAggregateColor,
    aggregateColor,
    poiDataset,
    workspace
  }: Props) => {
    const categoryOptions = poiDataset?.categories;

    React.useEffect(() => {
      // reset aggregate color to hide the aggregated dataset in the background
      setAggregateColor({key: "None", name: "None"})
    // eslint-disable-next-line
    }, [workspace]) // this is triggered during the embedding

    // React.useEffect(() => {
    //   if(poiDataset !== null && poiDataset.columns !== null){
    //     var catOpt = Object.keys(poiDataset.columns).map(key => {return {"key": key, "name": key};});
    //     setCategoryOptions({"attributes": catOpt});
    //   }
    // }, [poiDataset?.columns]);
    // console.log('AggregationTabPanel cimeBackgroundSelection:', cimeBackgroundSelection)


    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          {/* {
            categoryOptions != null ?
                <SelectFeatureComponent label={"color"} default_val={aggregateColor} categoryOptions={categoryOptions} onChange={(newValue) => {
                    var attribute = null
                    if (newValue && newValue != "") {
                        attribute = categoryOptions.attributes.filter(a => a.key == newValue)[0]
                    }
                    setAggregateColor(attribute)
                }}></SelectFeatureComponent>
                :
                <div></div>
        } */}
          {
            //TODO: check, if it makes sense to also include categorical values, or if it is ok to only use numerical values (like for "size")
            categoryOptions != null &&
            CategoryOptionsAPI.hasCategory(categoryOptions, "size") ? (
              <SelectFeatureComponent
                label={"color"}
                default_val={aggregateColor}
                categoryOptions={CategoryOptionsAPI.getCategory(
                  categoryOptions,
                  "size"
                )}
                onChange={(newValue) => {
                  var attribute = null;
                  if (newValue && newValue !== "") {
                    attribute = CategoryOptionsAPI.getCategory(
                      categoryOptions,
                      "size"
                    ).attributes.filter((a) => a.key === newValue)[0];
                  }
                  if (attribute === null || attribute === undefined) {
                    attribute = { key: "None", name: "None" };
                  }
                  setAggregateColor(attribute);

              }}></SelectFeatureComponent>)
              :
              <div></div>
          }
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          {
            <TextField
              id="knn-textfield"
              label="k-nearest neighbors"
              variant="outlined"
              defaultValue={50}
              fullWidth
              size="small"
              type="number"
            />
          }
        </Box>
      </div>
    );
  }
);

