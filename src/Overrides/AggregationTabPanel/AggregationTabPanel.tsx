import { TextField, Box } from "@mui/material";
import { CategoryOptionsAPI, SelectFeatureComponent, useCancellablePromise } from "projection-space-explorer";
// import { CategoryOptionsAPI, SelectFeatureComponent, useCancellablePromise } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import { setDatasetAction } from "projection-space-explorer/dist/components/Ducks/DatasetDuck"
import { setCimeBackgroundSelection } from "projection-space-explorer/dist/components/Ducks/CimeBackgroundSelectionDuck"
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { downloadImpl } from '../../Utility/Utils'
import React from "react";


const mapStateToProps = (state: AppState) => ({
  aggregateColor: state.aggregateColor,
  cimeBackgroundSelection: state.cimeBackgroundSelection,
  poiDataset: state.dataset,
  workspace: state.projections.workspace
});

const mapDispatchToProps = (dispatch) => ({
  setDataset: dataset => dispatch(setDatasetAction(dataset)),
  setAggregateColor: aggregateColor => dispatch(setAggregateColor(aggregateColor)),
  setCimeBackgroundSelection: (coords) => dispatch(setCimeBackgroundSelection(coords))

});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};

function hello(){
  console.log('debug button hello world')
}

export const AggregationTabPanel = connector(
  ({
    setAggregateColor,
    aggregateColor,
    poiDataset,
    cimeBackgroundSelection,
    setCimeBackgroundSelection,
    workspace
  }: Props) => {
    const { cancellablePromise, cancelPromises } = useCancellablePromise(); //TODO: cancelPromises --> use this to cancel promises on demand
    console.log("AgTabP.tsx poiDataset", poiDataset);
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

    handleBackgroundSelectionUpdate(cimeBackgroundSelection, setCimeBackgroundSelection, poiDataset?.info?.path);

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

/**
 * This is merely a helper function to decompose the code into smaller individual segments.
 * It is used to handle the changing background selection prop, which might, e.g., be triggered when a user clicks the k-nearest neighbor option in the context menu.
 * Specifically, it checks whether the parameters are correct, and if so, sends a query to the db to fetch the k-nearest entires to the click, with k being defined by a textfield.
 * The response triggers the download of a csv with these entries.
 * Afterwards, the background selection is reset, to make sure other prop updates do not trigger this db query and download.
 * @param {any} cimeBackgroundSelection - The prop form the redux store that holds x and y coordinates of the clicks
 * @param {any} setCimeBackgroundSelection - The dispatch function used to reset the prop after handling background selection
 * @returns {void} - no return value
 */
function handleBackgroundSelectionUpdate(cimeBackgroundSelection: any, setCimeBackgroundSelection: (coords: any) => any, filename: string) {
  // if input for checking k-nearest neighbors (x,y coordinates and k) are not undefined
  if (
    typeof cimeBackgroundSelection?.x !== "undefined" &&
    typeof cimeBackgroundSelection?.y !== "undefined" &&
    (document.getElementById("knn-textfield") as HTMLInputElement)?.value !==
      "undefined"
  ) {
    let k = +(document.getElementById("knn-textfield") as HTMLInputElement)
      ?.value;
    // if input k is neither integer nor below 1
    if (k < 1 || k % 1 !== 0) {
      // warn user
      alert("Invalid input for k-nearest neighbors.");
    } else {
      // otherwise send request to db and download response in browser
      ReactionCIMEBackendFromEnv.getkNearestData(
        filename,
        cimeBackgroundSelection?.x,
        cimeBackgroundSelection?.y,
        (document.getElementById("knn-textfield") as HTMLInputElement)?.value
      ).then((response) => {
        if (typeof response !== "undefined") {
          console.log(`response`, response);
          downloadImpl(
            JSON.stringify(response, null, 1),
            "k_nearest_data.csv",
            "text/csv"
          );
        }
      });
    }
    // reset x/y coordinates for safety, such that no other prop update will trigger this download
    setCimeBackgroundSelection(null);
  }
}

