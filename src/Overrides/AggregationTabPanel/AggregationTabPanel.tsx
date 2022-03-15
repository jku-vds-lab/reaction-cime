import { TextField, Box } from "@mui/material";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import React from "react";
import { StepSlider } from "./StepSlider";
import { ColorMapLegend } from "./ColorMapLegend";
import "./AggregationTabPanel.scss";
import { AdvancedAggregationSettings } from "./AdvancedAggregationSettings";
import { SelectFeatureComponent } from "projection-space-explorer";


const mapStateToProps = (state: AppState) => ({
  poiDataset: state.dataset,
  workspace: state.projections.workspace,
});

const mapDispatchToProps = (dispatch) => ({
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {};


export const AggregationTabPanel = connector(({poiDataset, workspace}: Props) => {
    
    // const categoryOptions = poiDataset?.categories;
    const [categoryOptions, setCategoryOptions] = React.useState(null)
    const [selectColumns, setSelectColumns] = React.useState(null)
    const [selectAttribute, setSelectAttribute] = React.useState({key:"None", name:"None"})
    const [selectAttributeInfo, setSelectAttributeInfo] = React.useState(null)
    const [columnInfo, setColumnInfo] = React.useState(null)
    
    React.useEffect(() => {
      // reset aggregate color to hide the aggregated dataset in the background
      setSelectAttribute({key:"None", name:"None"}) // reset selection to None
      setSelectAttributeInfo(null)
    // eslint-disable-next-line
    }, [workspace, poiDataset]) // this is triggered during the embedding


    React.useEffect(() => {
      setSelectAttribute({key:"None", name:"None"}) // reset selection to None
      setSelectAttributeInfo(null)

      let select_columns = {};
      let new_column_info = {...poiDataset?.columns};
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

                new_column_info[group_name] = {featureLabel: "Timeseries", distinct: undefined, featureType: undefined, isNumeric: undefined, metaInformation:undefined, project: undefined, range: undefined}
                
            }
            else if(col.metaInformation.real_column && col.isNumeric && !["x", "y"].includes(key)){ // we only want numeric features that exist in the real dataset and are not coordinates
              select_columns[key] = null
            }
        });

        setColumnInfo(new_column_info) // add columninfo for new timeseries columns
        setSelectColumns(select_columns)

        let cat_lst = Object.keys(select_columns).map(key => {
          return {key: key, name: key}
        });

        cat_lst.sort((valA, valB) => valA.name.localeCompare(valB.name))
        // cat_lst.splice(0,0,{key: "None", name: "None"})
        setCategoryOptions(cat_lst)

      }
    }, [poiDataset])


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
        //     <FormControl style={{ width: '100%' }}>
        //       <FormHelperText>Color by</FormHelperText>
        //       <Select 
        //       // fullWidth
        //         displayEmpty
        //         size='small'
        //         value={selectAttribute?.key}
        //         onChange={(event) => {
        //           let attribute = categoryOptions.filter(opt => opt.key === event.target.value)[0]
        //           if(attribute == null){
        //             attribute = {key: "None", name: "None"}
        //           }
                  
        //           setSelectAttribute(attribute)
        //           if(selectColumns != null && Object.keys(selectColumns).includes(attribute.key))
        //           setSelectAttributeInfo(selectColumns[attribute.key])
        //         }}
        //       >
        //       {categoryOptions.map(opt => { return <MenuItem key={opt.key} value={opt.key}>{opt.name}</MenuItem>})}
        //     </Select>
        // </FormControl>
        <SelectFeatureComponent
          column_info={columnInfo}
          label="color"
          default_val={selectAttribute}
          categoryOptions={{attributes: categoryOptions}}
          onChange={(newValue) => {
            let attribute = categoryOptions.filter(opt => opt.key === newValue)[0]
            if(attribute == null){
              attribute = {key: "None", name: "None"}
            }
            
            setSelectAttribute(attribute)
            if(selectColumns != null && Object.keys(selectColumns).includes(attribute.key))
            setSelectAttributeInfo(selectColumns[attribute.key])
          }}
      />
          }
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><StepSlider selectAttribute={selectAttribute} selectAttributeInfo={selectAttributeInfo}></StepSlider></Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><ColorMapLegend selectAttribute={selectAttribute}></ColorMapLegend></Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}><AdvancedAggregationSettings selectAttribute={selectAttribute} selectAttributeInfo={selectAttributeInfo}></AdvancedAggregationSettings></Box>
      </div>
    );
  }
);





