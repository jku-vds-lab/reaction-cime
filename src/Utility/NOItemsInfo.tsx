import { Typography } from "@mui/material";
import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { ReactionCIMEBackendFromEnv } from "../Backend/ReactionCIMEBackend";
import { AppState } from "../State/Store";

const mapStateToProps = (state: AppState) => ({
    dataset: state.dataset,
    globalLabels: state.globalLabels,
  });
  
  const mapDispatchToProps = (dispatch) => ({
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
    variant: "all" | "filterOutOfTotal"
  };
  
  export const NOItemsInfo = connector(({dataset, globalLabels, variant}: Props) => {
    const [totalDataPoints, setTotalDataPoints] = React.useState(-1)

    React.useEffect(()=> {
      if(dataset != null){
        ReactionCIMEBackendFromEnv.loadNODatapoints(dataset.info.path).then((res) => {
          setTotalDataPoints(res.no_datapoints)
        })
      }
    }, [dataset])

    let text = "";
    switch(variant){
      case "all":
        text = `Showing all ${totalDataPoints} ${globalLabels.itemLabelPlural} as aggregation`
        break;
      case "filterOutOfTotal":
        text = `Currently showing ${dataset?.vectors?.length} out of ${totalDataPoints} ${globalLabels.itemLabelPlural} (~${Math.round(dataset?.vectors?.length/totalDataPoints*10000)/100}% of the dataset)`
        break;

    }

    return totalDataPoints >= 0 && <Typography color="textSecondary" variant="body2">{text}</Typography>
  })