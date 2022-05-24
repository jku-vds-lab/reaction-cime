import React from "react";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { Checkbox, FormControlLabel, Grid, Radio, FormControl, InputLabel, Select, MenuItem, Typography } from "@mui/material";
import { AggregationMethod, setAggregationMethod, setDeriveRange, setUncertaintyRange, setValueRange, setVariableIndex, toggleDeriveRange } from "../../State/AggregateSettingsDuck";
import { MinMaxNumberInput } from "../../Utility/MinMaxNumberInput";
import { Box } from "@mui/system";

const mapStateToProps = (state: AppState) => ({
    aggregateSettings: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings,
    selectAttribute: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttribute,
    selectAttributeInfo: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttributeInfo,
});

const mapDispatchToProps = (dispatch) => ({
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
};
  

  
export const AdvancedAggregationSettings = connector(({ selectAttribute, setValueRange, setUncertaintyRange, toggleDeriveRange, 
            selectAttributeInfo, setVariableIndex, setAggregationMethod, aggregateSettings}: Props) => {
    if(selectAttribute == null || selectAttribute === "None"){
        return null;
    }
    
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
        <Box paddingTop={1}>
            <FormControlLabel control={<Checkbox checked={aggregateSettings.advancedSettings.deriveRange} onChange={()=>{toggleDeriveRange()}} title={"Derive Range from Data"}></Checkbox>} label="Derive Range from Data" />
            
            {(aggregateSettings.advancedSettings.valueRange != null && !aggregateSettings.advancedSettings.deriveRange) && 
                <MinMaxNumberInput title={`Customize Range for ` + aggregateSettings.colormapSettings.aggregateColor.value_col} target={"value"} range={aggregateSettings.advancedSettings.valueRange} setRange={setValueRange}></MinMaxNumberInput>
            }
            {(aggregateSettings.advancedSettings.uncertaintyRange != null && !aggregateSettings.advancedSettings.deriveRange) && 
                <MinMaxNumberInput title={`Customize Range for ` + aggregateSettings.colormapSettings.aggregateColor.uncertainty_col} target={"uncertainty"} range={aggregateSettings.advancedSettings.uncertaintyRange} setRange={setUncertaintyRange}></MinMaxNumberInput>
            }
        </Box>
        

        {selectAttributeInfo != null ? 
            <Box paddingTop={1} sx={{ flexGrow: 1 }}>
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
                                    setVariableIndex({valueVariableIndex: index, uncertaintyVariableIndex: aggregateSettings.advancedSettings.variableIndex.uncertaintyVariableIndex})
                                }
                            }} checked={index === aggregateSettings.advancedSettings.variableIndex.valueVariableIndex} /></Grid>
                        <Grid item xs={1}><Radio onChange={
                            (event) => {
                                if(event.target.checked){
                                    setVariableIndex({valueVariableIndex: aggregateSettings.advancedSettings.variableIndex.valueVariableIndex, uncertaintyVariableIndex: index})
                                }
                            }} checked={index === aggregateSettings.advancedSettings.variableIndex.uncertaintyVariableIndex} /></Grid>
                    </Grid>
                ))}
                <Grid container columns={{ xs: 3 }}>
                    <Grid item xs={1}></Grid>
                    <Grid item xs={1}>
                        <FormControl>
                            <InputLabel id="selectValueAggregation">Agg</InputLabel>
                            <Select
                                labelId="selectValueAggregation"
                                value={aggregateSettings.advancedSettings.aggregationMethod.valueAggregationMethod}
                                label="Value Aggregation"
                                onChange={(event)=>{setAggregationMethod({valueAggregationMethod: event.target.value, uncertaintyAggregationMethod: aggregateSettings.advancedSettings.aggregationMethod.uncertaintyAggregationMethod})}}
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
                                value={aggregateSettings.advancedSettings.aggregationMethod.uncertaintyAggregationMethod}
                                label="Uncertainty Aggregation"
                                onChange={(event)=>{setAggregationMethod({valueAggregationMethod: aggregateSettings.advancedSettings.aggregationMethod.valueAggregationMethod, uncertaintyAggregationMethod: event.target.value})}}
                            >
                                {Object.values(AggregationMethod).map((value) => <MenuItem value={value} key={value}>{value}</MenuItem>)}
                            </Select>
                        </FormControl>
                    </Grid>
                </Grid>
                </>
            </Box>:
            <Box paddingTop={1}>
                <Typography variant="subtitle2" gutterBottom>Select value aggregation function</Typography>
                <FormControl>
                    <InputLabel id="selectValueAggregation">Agg</InputLabel>
                    <Select
                        labelId="selectValueAggregation"
                        value={aggregateSettings.advancedSettings.aggregationMethod.valueAggregationMethod}
                        label="Value Aggregation"
                        onChange={(event)=>{setAggregationMethod({valueAggregationMethod: event.target.value, uncertaintyAggregationMethod: aggregateSettings.advancedSettings.aggregationMethod.uncertaintyAggregationMethod})}}
                    >
                        {Object.values(AggregationMethod).map((value) => <MenuItem value={value} key={value}>{value}</MenuItem>)}
                    </Select>
                </FormControl>
            </Box>
        }
    </>
      
})
  
  