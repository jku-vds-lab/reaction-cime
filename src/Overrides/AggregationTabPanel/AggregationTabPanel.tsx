import { TextField, Box, Select, MenuItem, FormControl, FormHelperText } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import React from "react";
import { StepSlider } from "./StepSlider";
import { ColorMapLegend } from "./ColorMapLegend";


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

        let cat_lst = Object.keys(select_columns).map(key => {
          return {key: key, name: key, col_info: select_columns[key]}
        });

        cat_lst.sort((valA, valB) => valA.name.localeCompare(valB.name))
        cat_lst.splice(0,0,{key: "None", name: "None", col_info:null})
        setCategoryOptions(cat_lst)

      }
    }, [poiDataset])



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
        
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><StepSlider selectAttribute={selectAttribute}></StepSlider></Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><ColorMapLegend selectAttribute={selectAttribute}></ColorMapLegend></Box>
      </div>
    );
  }
);





