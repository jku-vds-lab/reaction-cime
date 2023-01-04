import { Box } from '@mui/material';
import { connect, ConnectedProps } from 'react-redux';
import React from 'react';
import { SelectFeatureComponent } from 'projection-space-explorer';
import { AppState } from '../../State/Store';
import { StepSlider } from './StepSlider';
import { ColorMapLegend } from './ColorMapLegend';
import './AggregationTabPanel.scss';
import { AdvancedAggregationSettings } from './AdvancedAggregationSettings';
import { AggregateActions } from '../../State/AggregateSettingsDuck';
import { NOItemsInfo } from '../../Utility/NOItemsInfo';

const mapStateToProps = (state: AppState) => ({
  poiDataset: state.dataset,
  selectAttribute: state.multiples.multiples.entities[state.multiples.active]?.attributes.aggregateSettings?.selectAttribute,
  globalLabels: state.globalLabels,
});

const mapDispatchToProps = (dispatch) => ({
  setSelectAttribute: (name: string, info) => dispatch(AggregateActions.setSelectAttribute({ attribute_name: name, attribute_info: info })),
  setValueRange: (range) => dispatch(AggregateActions.setValueRange(range)),
  setUncertaintyRange: (range) => dispatch(AggregateActions.setUncertaintyRange(range)),
  setDeriveRange: (value) => dispatch(AggregateActions.setDeriveRange(value)),
});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux;

export const AggregationTabPanel = connector(
  ({ poiDataset, selectAttribute, setSelectAttribute, setValueRange, setUncertaintyRange, setDeriveRange }: Props) => {
    // const categoryOptions = poiDataset?.categories;
    const [categoryOptions, setCategoryOptions] = React.useState(null);
    const [selectColumns, setSelectColumns] = React.useState(null);
    // const [selectAttribute, setSelectAttribute] = React.useState(null)
    // const [selectAttributeInfo, setSelectAttributeInfo] = React.useState(null)
    const [columnInfo, setColumnInfo] = React.useState(null);

    React.useEffect(() => {
      setValueRange(null);
      setUncertaintyRange(null);
      setDeriveRange(true);

      // eslint-disable-next-line
    }, [selectAttribute]);

    React.useEffect(() => {
      // setSelectAttribute(null) // reset selection to None
      // setSelectAttributeInfo(null)

      const newSelectColumns = {};
      const newColumnInfo = { ...poiDataset?.columns };
      if (poiDataset != null && poiDataset.columns != null) {
        Object.keys(poiDataset.columns).forEach((key) => {
          const col = poiDataset.columns[key];
          if (Object.keys(col.metaInformation).length > 0 && col.metaInformation.timeSeriesGroup) {
            const split = col.metaInformation.timeSeriesGroup.split(':');

            let groupName;
            let varName;
            if (split.length <= 1) {
              // if the string is separated with a colon, only the first part of the string is considered as the group. the second part of the string determines a sub value of this group
              groupName = col.metaInformation.timeSeriesGroup;
              varName = col.metaInformation.timeSeriesGroup;
            } else {
              groupName = split[0];
              varName = split[1];
            }
            // TODO: timestep information -> what if it is not available? extract from column name?
            // TODO: might want to have mean/variance values for non-timestap values bzw. single values
            if (Object.keys(newSelectColumns).includes(groupName)) {
              if (Object.keys(newSelectColumns[groupName]).includes(varName)) {
                newSelectColumns[groupName][varName].temporal_columns[col.metaInformation.timestep] = key;
              } else {
                const tempCols = {};
                tempCols[col.metaInformation.timestep] = key;
                newSelectColumns[groupName][varName] = { temporal_columns: tempCols };
              }
            } else {
              const tempCols = {};
              tempCols[col.metaInformation.timestep] = key;
              newSelectColumns[groupName] = {};
              newSelectColumns[groupName][varName] = { temporal_columns: tempCols };
            }

            newColumnInfo[groupName] = {
              featureLabel: 'Timeseries',
              distinct: undefined,
              featureType: undefined,
              isNumeric: undefined,
              metaInformation: undefined,
              project: undefined,
              range: undefined,
            };
          } else if (col.metaInformation.real_column && col.isNumeric && !['x', 'y'].includes(key)) {
            // we only want numeric features that exist in the real dataset and are not coordinates
            newSelectColumns[key] = null;
          }
        });

        setColumnInfo(newColumnInfo); // add columninfo for new timeseries columns
        setSelectColumns(newSelectColumns);

        // let cat_lst = Object.keys(select_columns).map(key => {
        //   return {key: key, name: key}
        // });
        const catLst = Object.keys(newSelectColumns);
        catLst.sort();
        // cat_lst.sort((valA, valB) => valA.name.localeCompare(valB.name))
        // cat_lst.splice(0,0,{key: "None", name: "None"})
        setCategoryOptions(catLst);
      }
    }, [poiDataset]);

    return (
      <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <NOItemsInfo variant="all" />
        </Box>
        {/* <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
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
        </Box> */}

        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          {
            // TODO: check, if it makes sense to also include categorical values, or if it is ok to only use numerical values (like for "size")
            categoryOptions != null && (
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
                categoryOptions={categoryOptions}
                onChange={(newValue) => {
                  let attribute = newValue; // categoryOptions.filter(opt => opt.key === newValue)[0]
                  if (attribute == null) {
                    attribute = 'None'; // {key: "None", name: "None"}
                  }
                  // setSelectAttribute(attribute)
                  if (selectColumns != null && Object.keys(selectColumns).includes(attribute)) {
                    setSelectAttribute(attribute, selectColumns[attribute]);
                    // setSelectAttributeInfo(selectColumns[attribute])
                  } else {
                    setSelectAttribute(attribute, null);
                  }
                }}
              />
            )
          }
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <StepSlider />
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <ColorMapLegend />
        </Box>
        <Box paddingLeft={2} paddingTop={1} paddingRight={2}>
          <AdvancedAggregationSettings />
        </Box>
      </div>
    );
  },
);
