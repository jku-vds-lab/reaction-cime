import { MenuItem } from "@mui/material"
import { CameraTransformations, IProjection, TypedObject } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { AppState } from "../../State/Store";

const mapStateToProps = (state: AppState) => ({
    dataset: state.dataset,
    multiples: state.multiples,
    globalLabels: state.globalLabels,
    triggerDatasetUpdate: state.handleDataset?.triggerUpdate,
    state: state
  });
  
  const mapDispatchToProps = (dispatch) => ({
  });
  
  const connector = connect(mapStateToProps, mapDispatchToProps);
  
  type PropsFromRedux = ConnectedProps<typeof connector>;
  
  type Props = PropsFromRedux & {
    handleClose: () => void, 
    pos_x: number, 
    pos_y: number,
    menuTarget: TypedObject,
  };
  
export const AddRegionExceptionMenuItem = connector(({ handleClose, pos_x, pos_y, dataset, multiples, globalLabels, triggerDatasetUpdate, state }: Props) => {
    const workspace = multiples.multiples.entities[multiples.active]?.attributes?.workspace as IProjection;
    const xChannel = workspace?.xChannel == null ? "x" : workspace?.xChannel;
    const yChannel = workspace?.yChannel == null ? "y" : workspace?.yChannel;
    
    const viewTransform = multiples.multiples.entities[multiples.active]?.attributes?.viewTransform;
    
    return <MenuItem
      onClick={() => {
        const coords = CameraTransformations.screenToWorld({x: pos_x, y: pos_y}, viewTransform);

        const deltaX = dataset.columns[xChannel].range.max - dataset.columns[xChannel].range.min;
        const deltaY = dataset.columns[yChannel].range.max - dataset.columns[yChannel].range.min;
        const radius = Math.min(deltaX, deltaY)*0.1;

        ReactionCIMEBackendFromEnv.addPOIExceptions(dataset?.info.path, [{x_col: xChannel, y_col: yChannel, x_coord: coords.x, y_coord: coords.y, radius: radius}]).then((res) => {
            if(res.msg !== "ok"){
                alert(res.msg)
            }
            if(triggerDatasetUpdate != null){
                triggerDatasetUpdate({
                    display: dataset.info.path,
                    path: dataset.info.path,
                    type: dataset.info.type,
                    uploaded: true
                }, state)
            }
            
        })
        
        handleClose()
      }}
    >
      {`Show ${globalLabels.itemLabelPlural} in this region`}
    </MenuItem>
  })