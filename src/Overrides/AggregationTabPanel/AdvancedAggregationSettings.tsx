import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { Button, Checkbox, FormControlLabel, TextField, Grid, Radio, FormControl, InputLabel, Select, MenuItem } from "@mui/material";
import { AggregationMethod, setAggregationMethod, setDeriveRange, setSampleSize, setUncertaintyRange, setValueRange, setVariableIndex, toggleDeriveRange } from "../../State/AggregateSettingsDuck";
import { MinMaxNumberInput } from "../../Utility/MinMaxNumberInput";
import { Box } from "@mui/system";

const mapStateToProps = (state: AppState) => ({
    sampleSize: state.aggregateSettings?.sampleSize,
    aggregateColor: state.aggregateSettings?.aggregateColor,
    valueRange: state.aggregateSettings?.valueRange,
    uncertaintyRange: state.aggregateSettings?.uncertaintyRange,
    deriveRange: state.aggregateSettings?.deriveRange,
    variableIndex: state.aggregateSettings?.variableIndex,
    aggregationMethod: state.aggregateSettings?.aggregationMethod
});

const mapDispatchToProps = (dispatch) => ({
    setSampleSize: value => dispatch(setSampleSize(value)),
    setValueRange: (range) => dispatch(setValueRange(range)),
    setUncertaintyRange: (range) => dispatch(setUncertaintyRange(range)),
    toggleDeriveRange: () => dispatch(toggleDeriveRange()),
    setDeriveRange: (value) => dispatch(setDeriveRange(value)),
    setVariableIndex: (value) => dispatch(setVariableIndex(value)),
    setAggregationMethod: (value) => dispatch(setAggregationMethod(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    selectAttribute,
    selectAttributeInfo
};
  

  
export const AdvancedAggregationSettings = connector(({sampleSize, setSampleSize, selectAttribute, aggregateColor, 
            valueRange, uncertaintyRange, setValueRange, setUncertaintyRange, deriveRange, toggleDeriveRange, 
            setDeriveRange, selectAttributeInfo, setVariableIndex, setAggregationMethod, variableIndex, aggregationMethod}: Props) => {
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


    React.useEffect(()=>{
        if(selectAttributeInfo){
            // initialize valueIndices when attributionInfo changes
            const variables_array = Object.keys(selectAttributeInfo);
            if(variables_array.length === 1){
                setVariableIndex({valueVariableIndex: 0, uncertaintyVariableIndex: 0})
            }else if(variables_array.length > 1){
                setVariableIndex({valueVariableIndex: 0, uncertaintyVariableIndex: 1})
            }
        }
    // eslint-disable-next-line
    }, [selectAttributeInfo])


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

        {selectAttributeInfo != null ? 
            <Box sx={{ flexGrow: 1 }}>
                <>
                    <Grid container columns={{ xs: 3 }}>
                        <Grid item xs={1}></Grid>
                        <Grid item xs={1}>value</Grid>
                        <Grid item xs={1}>uncertainty</Grid>
                    </Grid>
                {Object.keys(selectAttributeInfo).map((value, index) => (
                    <Grid container columns={{ xs: 3 }} key={value}>
                        <Grid item xs={1}>{value}</Grid>
                        <Grid item xs={1}><Radio onChange={
                            (event) => {
                                if(event.target.checked){
                                    setVariableIndex({valueVariableIndex: index, uncertaintyVariableIndex: variableIndex.uncertaintyVariableIndex})
                                }
                            }} checked={index === variableIndex.valueVariableIndex} /></Grid>
                        <Grid item xs={1}><Radio onChange={
                            (event) => {
                                if(event.target.checked){
                                    setVariableIndex({valueVariableIndex: variableIndex.valueVariableIndex, uncertaintyVariableIndex: index})
                                }
                            }} checked={index === variableIndex.uncertaintyVariableIndex} /></Grid>
                    </Grid>
                ))}
                <Grid container columns={{ xs: 3 }}>
                    <Grid item xs={1}></Grid>
                    <Grid item xs={1}>
                        <FormControl>
                            <InputLabel id="selectValueAggregation">Agg</InputLabel>
                            <Select
                                labelId="selectValueAggregation"
                                value={aggregationMethod.valueAggregationMethod}
                                label="Value Aggregation"
                                onChange={(event)=>{setAggregationMethod({valueAggregationMethod: event.target.value, uncertaintyAggregationMethod: aggregationMethod.uncertaintyAggregationMethod})}}
                            >
                                {Object.values(AggregationMethod).map((value) => <MenuItem value={value} key={value}>{value}</MenuItem>)}
                            </Select>
                        </FormControl>
                    </Grid>
                    <Grid item xs={1}>
                        
                    <FormControl>
                            <InputLabel id="selectUncertaintyAggregation">Agg</InputLabel>
                            <Select
                                labelId="selectUncertaintyAggregation"
                                value={aggregationMethod.uncertaintyAggregationMethod}
                                label="Uncertainty Aggregation"
                                onChange={(event)=>{setAggregationMethod({valueAggregationMethod: aggregationMethod.valueAggregationMethod, uncertaintyAggregationMethod: event.target.value})}}
                            >
                                {Object.values(AggregationMethod).map((value) => <MenuItem value={value} key={value}>{value}</MenuItem>)}
                            </Select>
                        </FormControl>
                    </Grid>
                </Grid>
                </>
            </Box>:
            <Box>
                <FormControl>
                    <InputLabel id="selectValueAggregation">Agg</InputLabel>
                    <Select
                        labelId="selectValueAggregation"
                        value={aggregationMethod.valueAggregationMethod}
                        label="Value Aggregation"
                        onChange={(event)=>{setAggregationMethod({valueAggregationMethod: event.target.value, uncertaintyAggregationMethod: aggregationMethod.uncertaintyAggregationMethod})}}
                    >
                        {Object.values(AggregationMethod).map((value) => <MenuItem value={value} key={value}>{value}</MenuItem>)}
                    </Select>
                </FormControl>
            </Box>
        }
    </>
      
})
  
  