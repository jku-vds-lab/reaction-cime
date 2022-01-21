import { TextField, Box, Select, MenuItem, FormControl, FormHelperText, Typography, Slider } from "@mui/material";
// import { CategoryOptionsAPI, SelectFeatureComponent, useCancellablePromise } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { setAggregateColor } from "../../State/AggregateColorDuck";
import React from "react";
import * as d3 from 'd3v5';


const mapStateToProps = (state: AppState) => ({
  poiDataset: state.dataset,
  legend: state.aggregateSettings?.legend
});

const mapDispatchToProps = (dispatch) => ({
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};


export const AggregationTabPanel = connector(
  ({
    poiDataset,
    legend
  }: Props) => {
    
    // const categoryOptions = poiDataset?.categories;
    const [categoryOptions, setCategoryOptions] = React.useState(null)
    const [selectAttribute, setSelectAttribute] = React.useState({key:"None", name:"None", col_info:null})


    React.useEffect(() => {

      let select_columns = {};
      if(poiDataset != null && poiDataset.columns != null){
        Object.keys(poiDataset.columns).forEach(key => {
          let col = poiDataset.columns[key];
            if(Object.keys(col.metaInformation).length > 0 && col.metaInformation.timeSeriesGroup){
                const split = col.metaInformation.timeSeriesGroup.split(":"); 
                
                let group_name, var_name;
                if(split.length <=1){ // if the string is separated with a colon, only the first part of the string is considered as the group. the second part of the string determines a sub value of this group
                  group_name = col.metaInformation.timeSeriesGroup;
                  var_name = col.metaInformation.timeSeriesGroup;
                }else{
                  group_name = split[0];
                  var_name = split[1];
                }
                // TODO: timestep information -> what if it is not available? extract from column name?
                // TODO: might want to have mean/variance values for non-timestap values bzw. single values
                if(Object.keys(select_columns).includes(group_name)){
                  if(Object.keys(select_columns[group_name]).includes(var_name)){
                    select_columns[group_name][var_name].temporal_columns[col.metaInformation.timestep] = key
                  }else{
                    let temp_cols = {};
                    temp_cols[col.metaInformation.timestep] = key
                    select_columns[group_name][var_name] = {"temporal_columns": temp_cols} 
                  }
                }else{
                  let temp_cols = {};
                  temp_cols[col.metaInformation.timestep] = key
                  select_columns[group_name] = {}
                  select_columns[group_name][var_name] = {"temporal_columns": temp_cols} 
                }
                
            }
            // else{
            //     row[i] = element[i];
            //     columns[i] = col;
            // }
        });

        console.log(select_columns)
        let cat_lst = Object.keys(select_columns).map(key => {
          return {key: key, name: key, col_info: select_columns[key]}
        });

        cat_lst.sort((valA, valB) => valA.name.localeCompare(valB.name))
        cat_lst.splice(0,0,{key: "None", name: "None", col_info:null})
        setCategoryOptions(cat_lst)

      }
    }, [poiDataset])

    // add colormap legend
    const gRef = React.useRef();
    React.useEffect(() => {
      if(legend != null){
        const gElement = d3.select(gRef.current)
        gElement.html(""); // clear g element
        gElement.call(legend); // draw color legend
      }
      // const svgElement = d3.select(svgRef.current)
      // var vDom = [0,100];
      // var uDom = [0,1];

      // var quantization = vsup.quantization().branching(3).layers(3).valueDomain(vDom).uncertaintyDomain(uDom);
      // var scale = vsup.scale().quantize(quantization).range(d3.interpolateYlGnBu);
      // var legend = vsup.legend.arcmapLegend(scale);
      // legend
      // //   .scale(scale)
      //   .size(180)
      //   .x(20)
      //   .y(50)
      //   .vtitle("Mean Prediction")
      //   .utitle("Variance in Prediction");

      // svgElement.append("g").call(legend)
    }, [legend])


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

        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          {
            //TODO: check, if it makes sense to also include categorical values, or if it is ok to only use numerical values (like for "size")
            categoryOptions != null && 
            <FormControl style={{ width: '100%' }}>
              <FormHelperText>Color by</FormHelperText>
              <Select 
              // fullWidth
                displayEmpty
                size='small'
                value={selectAttribute?.key}
                onChange={(event) => {
                  let attribute = categoryOptions.filter(opt => opt.key === event.target.value)[0]
                  if(attribute == null){
                    attribute = {key: "None", name: "None", col_info: null}
                  }
                  
                  setSelectAttribute(attribute)
                }}
              >
              {categoryOptions.map(opt => { return <MenuItem key={opt.key} value={opt.key}>{opt.name}</MenuItem>})}
            </Select>
        </FormControl>
          }
        </Box>
        {/* <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          {
            //TODO: check, if it makes sense to also include categorical values, or if it is ok to only use numerical values (like for "size")
            categoryOptions != null &&
            CategoryOptionsAPI.hasCategory(categoryOptions, "size") ? (
              <SelectFeatureComponent
                column_info={poiDataset?.columns}
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
        </Box> */}
        
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><svg style={{width:"100%", height:"300px"}}><g ref={gRef}></g></svg></Box>
        <Box><StepSlider selectAttribute={selectAttribute}></StepSlider></Box>
      </div>
    );
  }
);

const mapStateToPropsSlider = (state: AppState) => ({
  workspace: state.projections.workspace,
  poiDataset: state.dataset
});

const mapDispatchToPropsSlider = (dispatch) => ({
  setAggregateColor: values => dispatch(setAggregateColor(values)),
});

const sliderconnector = connect(mapStateToPropsSlider, mapDispatchToPropsSlider);

type SliderPropsFromRedux = ConnectedProps<typeof sliderconnector>;

type SliderProps = SliderPropsFromRedux & {
  selectAttribute: {key:string, name:string, col_info:any}
};



const StepSlider = sliderconnector(({selectAttribute, setAggregateColor, workspace, poiDataset}: SliderProps) => {
  if(selectAttribute == null || selectAttribute.key === "None" || selectAttribute.key == null){
    return null;
  }
  
  React.useEffect(() => {
    // reset aggregate color to hide the aggregated dataset in the background
    setAggregateColor(null)
  // eslint-disable-next-line
  }, [workspace, poiDataset]) // this is triggered during the embedding

  
  const stepChanged = (newVal) => {
    setCurStep(newVal)
  }
  
  const [curStep, setCurStep] = React.useState(0);
  const [marks, setMarks] = React.useState([]);
  

  React.useEffect(() => {
    let m;
    if(selectAttribute.col_info){
      const variables = Object.keys(selectAttribute.col_info);
      let step_arr = Object.keys(selectAttribute.col_info[variables[0]].temporal_columns)
      step_arr.sort()

      m = step_arr.map((step, index) => {return {value: step, label: step}})
      console.log(m)
    }
    setMarks(m)
    setCurStep(m[m.length-1].value)
  }, [selectAttribute])

  React.useEffect(() => {
    const timestep = curStep;
    if(selectAttribute.col_info){
      const variables_array = Object.keys(selectAttribute.col_info);
      // TODO: users should be able to select which column is value and which is uncertainty
      let value_col = selectAttribute.col_info[variables_array[0]]["temporal_columns"][timestep]
      let uncertainty_col = null;
      if(variables_array.length >= 2){
        uncertainty_col = selectAttribute.col_info[variables_array[1]]["temporal_columns"][timestep]
      }
      setAggregateColor({"value_col": value_col, "uncertainty_col": uncertainty_col});
    }else{
      setAggregateColor({"value_col": selectAttribute.key, "uncertainty_col": null});
    }
    // eslint-disable-next-line
  }, [curStep])

  
  
  

  return marks && marks.length > 0 && <div style={{
      margin: '0px 16px',
      padding: '0px 8px'
  }}>
      <Typography id="range-slider" gutterBottom>
          Choose Step
      </Typography>
      <Slider
          min={0}
          max={marks.length-1}
          value={curStep}
          onChange={(_, newValue) => stepChanged(newValue)}
          step={1}
          marks={marks}
          valueLabelDisplay="auto"
      ></Slider>
  </div>
})