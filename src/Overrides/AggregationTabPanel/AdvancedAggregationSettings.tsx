import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { Button, TextField } from "@mui/material";
import { setSampleSize } from "../../State/AggregateSettingsDuck";

const mapStateToProps = (state: AppState) => ({
    sampleSize: state.aggregateSettings?.sampleSize,
});

const mapDispatchToProps = (dispatch) => ({
    setSampleSize: value => dispatch(setSampleSize(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
};
  

  
export const AdvancedAggregationSettings = connector(({sampleSize, setSampleSize}: Props) => {
    const [tempSampleSize, setTempSampleSize] = React.useState(sampleSize);

    React.useEffect(() => {
        setTempSampleSize(sampleSize)
    }, [sampleSize])
    
    const handleChange = (event) => {
        let val = event.target.value
        val = val > 300 ? 300 : val
        val = val < 5 ? 5 : val
        setTempSampleSize(val);
    };

    // TODO: should we give a user input for sample size? the benefits are not so big (would be a maximum of 300 due to browser memory issues, now it is set to 200)
    // if we decide to include it: remove "false" flag; also make sure that cache is cleared and that the background is updated when sample size is changed 
    return <>{false && <> 
        <TextField
          id="sample-size"
          
          label="Sample Size"
          type="number"
          value={tempSampleSize}
          InputLabelProps={{
            shrink: true,
          }}
          onChange={handleChange}
        />
        <Button onClick={() => {setSampleSize(tempSampleSize)}}>Apply Settings</Button>
    </>}</>
      
})
  
  