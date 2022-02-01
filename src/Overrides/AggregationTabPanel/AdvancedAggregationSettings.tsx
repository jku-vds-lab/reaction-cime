import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { Button, Checkbox, FormControlLabel, TextField } from "@mui/material";
import { setDeriveRange, setSampleSize, setUncertaintyRange, setValueRange, toggleDeriveRange } from "../../State/AggregateSettingsDuck";
import { MinMaxNumberInput } from "../../Utility/MinMaxNumberInput";

const mapStateToProps = (state: AppState) => ({
    sampleSize: state.aggregateSettings?.sampleSize,
    aggregateColor: state.aggregateSettings?.aggregateColor,
    valueRange: state.aggregateSettings?.valueRange,
    uncertaintyRange: state.aggregateSettings?.uncertaintyRange,
    deriveRange: state.aggregateSettings?.deriveRange,
});

const mapDispatchToProps = (dispatch) => ({
    setSampleSize: value => dispatch(setSampleSize(value)),
    setValueRange: (range) => dispatch(setValueRange(range)),
    setUncertaintyRange: (range) => dispatch(setUncertaintyRange(range)),
    toggleDeriveRange: () => dispatch(toggleDeriveRange()),
    setDeriveRange: (value) => dispatch(setDeriveRange(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    selectAttribute
};
  

  
export const AdvancedAggregationSettings = connector(({sampleSize, setSampleSize, selectAttribute, aggregateColor, valueRange, uncertaintyRange, setValueRange, setUncertaintyRange, deriveRange, toggleDeriveRange, setDeriveRange}: Props) => {
    if(selectAttribute == null || selectAttribute.key === "None" || selectAttribute.key == null){
        return null;
    }

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

    React.useEffect(() => {
        setValueRange(null)
        setUncertaintyRange(null)
        setDeriveRange(true)

        // eslint-disable-next-line
    }, [selectAttribute])



    // TODO: should we give a user input for sample size? the benefits are not so big (would be a maximum of 300 due to browser memory issues, now it is set to 200)
    // if we decide to include it: remove "false" flag; also make sure that cache is cleared and that the background is updated when sample size is changed 
    return <>
        {false && <> 
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
        </>}
        <FormControlLabel control={<Checkbox checked={deriveRange} onChange={()=>{toggleDeriveRange()}} title={"Derive Range from Data"}></Checkbox>} label="Derive Range from Data" />
        
        {(valueRange != null && !deriveRange) && 
            <MinMaxNumberInput title={`Customize Range for ` + aggregateColor.value_col} target={"value"} range={valueRange} setRange={setValueRange}></MinMaxNumberInput>
        }
        {(uncertaintyRange != null && !deriveRange) && 
            <MinMaxNumberInput title={`Customize Range for ` + aggregateColor.uncertainty_col} target={"uncertainty"} range={uncertaintyRange} setRange={setUncertaintyRange}></MinMaxNumberInput>
        }
    </>
      
})
  
  