import { Box } from "@mui/material";
import { CategoryOptionsAPI, SelectFeatureComponent, useCancellablePromise } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { setAggregateDatasetAction } from "../../State/AggregateDatasetDuck";
import { AppState } from "../../State/Store";
import { AggregateDataset } from "./AggregateDataset";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import dataset from "projection-space-explorer/dist/components/Ducks/DatasetDuck";

const mapStateToProps = (state: AppState) => ({
  aggregateColor: state.aggregateColor,
  poiDataset: state.dataset,
});

const mapDispatchToProps = (dispatch) => ({
  setAggregateDataset: dataset => dispatch(setAggregateDatasetAction(dataset)),
  setAggregateColor: aggregateColor => dispatch(setAggregateColor(aggregateColor)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};

export const AggregationTabPanel = connector(({setAggregateDataset, setAggregateColor, aggregateColor, poiDataset}: Props) => {
    const { cancellablePromise, cancelPromises } = useCancellablePromise(); //TODO: cancelPromises --> use this to cancel promises on demand
    console.log('AgTabP.tsx poiDataset', poiDataset)
    const categoryOptions = poiDataset?.categories;

    // let [categoryOptions, setCategoryOptions] = React.useState(null);

    // React.useEffect(() => {
    //   if(poiDataset !== null && poiDataset.columns !== null){
    //     var catOpt = Object.keys(poiDataset.columns).map(key => {return {"key": key, "name": key};});
    //     setCategoryOptions({"attributes": catOpt});
    //   }
    // }, [poiDataset?.columns]);
    
    return (
      <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
        {/* <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <Typography variant="subtitle2" gutterBottom>
            Upload Aggregated Dataset
          </Typography>
        </Box>
        <AggregatedDatasetDrop onChange={(dataset) => {setAggregateDataset(new AggregateDataset(dataset))}} /> */}
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

                  if(attribute.key === "None"){
                    setAggregateDataset(null)
                  } else {
                    let abort_controller = new AbortController(); //TODO: reiterate where AbortController needs to be instantiated --> can it be moved inside the loadAggCSV function?
                    ReactionCIMEBackendFromEnv.loadAggCSV((dataset) => {setAggregateDataset(new AggregateDataset(dataset))}, poiDataset.info.path, attribute.key, cancellablePromise, null, abort_controller)
                  }
              }}></SelectFeatureComponent>
              :
              <div></div>
        }
        </Box>
      </div>
    );
  }
);
