import { Box, Button } from "@mui/material";
import { Dataset } from "projection-space-explorer";
import React from "react";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { CategoryFilter } from "./CategoryFilter";
import { RangeFilter } from "./RangeFilter";
import FilterAltIcon from '@mui/icons-material/FilterAlt';



type Props = {
    dataset: Dataset,
    removeFilter: (col) => void,
    // constraintsMap: {}
    constraintCols: string[],
    constraints: {col:string, operator:string, val1:string, val2:string}[],
    triggerDatasetUpdate
};

export const FilterSettings = ({dataset, removeFilter, constraintCols, constraints, triggerDatasetUpdate}:Props) => {

    const [filterValues, setFilterValues] = React.useState({});

    React.useEffect(() => {
        // initialize filters
        if(constraints != null){
            let tempFilterValues = {...filterValues};
            for(let i in constraintCols){
                const col = constraintCols[i];
                const currentConstraints = constraints.filter((item) => item.col === col);

                let between = [undefined, undefined];
                let equals = [];
                if(currentConstraints.length > 0){
                    currentConstraints.forEach((item) => {
                        for(let i in constraints){
                            const con = constraints[i];
                            if(con.operator === "BETWEEN"){ // if there are several between operators, we choose the minimal minimum and the maximal maximum, because we can only handle one range
                                if(between[0] == null || between[1] == null){
                                    between = [parseFloat(con.val1), parseFloat(con.val2)]
                                }else{
                                    between[0] = Math.min(parseFloat(con.val1), between[0])
                                    between[1] = Math.max(parseFloat(con.val2), between[1])
                                }
                                
                            }else if(con.operator === "EQUALS"){
                                equals.push(con.val1)
                            }
                        }
                    })
                }
                
                if(dataset.columns[col].isNumeric){
                    tempFilterValues[col] = {isNum: true, val: between}
                }else{
                    tempFilterValues[col] = {isNum: false, val: equals}
                }
            }
            setFilterValues(tempFilterValues);
        }
        
    }, [constraints])

    React.useEffect(() => {
        if(constraintCols != null){
            let tempFilterValues = {...filterValues};
            for(let i in constraintCols){
                // add constraints, if not already included
                const col = constraintCols[i];
                if(!Object.keys(tempFilterValues).includes(col)){
                    if(dataset.columns[col].isNumeric){
                        tempFilterValues[col] = {isNum: true, val: [-Number.MAX_VALUE, Number.MAX_VALUE]}
                    }else{
                        tempFilterValues[col] = {isNum: false, val: []}
                    }
                }
            }
            for(let i in Object.keys(tempFilterValues)){
                // remove constraints
                const col = Object.keys(tempFilterValues)[i];
                if(!constraintCols.includes(col)){
                    delete tempFilterValues[col];
                }
            }
            setFilterValues(tempFilterValues);
        }
    }, [constraintCols])

    // React.useEffect(()=> {
    //     Object.keys(constraintsMap).map((key) => {
    //         const constraints = constraintsMap[key];
    //         let between;
    //         let equals = [];
    //         for(let i in constraints){
    //             const con = constraints[i];
    //             if(con.operator === "BETWEEN"){ // if there are several between operators, we choose the minimal minimum and the maximal maximum, because we can only handle one range
    //                 if(between == null){
    //                     between = [con.val1, con.val2]
    //                 }else{
    //                     between[0] = Math.min(con.val1, between[0])
    //                     between[1] = Math.max(con.val2, between[1])
    //                 }
                    
    //             }else if(con.operator === "EQUALS"){
    //                 equals.push(con.val1)
    //             }
    //         }
    //         let tempFilterValues = {...filterValues};
    //         if(between != null){
    //             tempFilterValues[key] = {isNum: true, val: between}
    //         }else{
    //             tempFilterValues[key] = {isNum: false, val: equals}
    //         }
    //         setFilterValues(tempFilterValues);
    //         return null;
    //     })

    // }, [constraintsMap])

    return <div><Box paddingTop={2} paddingRight={2}>
        {Object.keys(filterValues).map((key) => {
            const value = filterValues[key]

            if(value.isNum){
                return <RangeFilter key={key} dataset={dataset} col={key} value={value.val} setValue={(newValue) => {
                    let tempFilterValues = {...filterValues};
                    tempFilterValues[key].val = newValue
                    setFilterValues(tempFilterValues);
                }} remove={removeFilter}></RangeFilter>
            }else{
                return <CategoryFilter key={key} dataset={dataset} col={key} value={value.val} setValue={(newValue) => {
                    let tempFilterValues = {...filterValues};
                    tempFilterValues[key].val = newValue
                    setFilterValues(tempFilterValues);
                }} remove={removeFilter}></CategoryFilter>
            }
        })}
        
        </Box>
        
        <Box paddingLeft={2} paddingTop={1}>
            <Button
                fullWidth
                variant="outlined"
                onClick={() => {
                    updateBackendConstraints(filterValues, dataset, triggerDatasetUpdate)
                }}
            >
            <FilterAltIcon />
            &nbsp;Apply Filter
            </Button>
        </Box>
        </div>
}


const updateBackendConstraints = (dimensions, dataset, triggerDatasetUpdate) => {
    const constraint_dimensions = dimensions;
    let all_constraints = []
    for(const i in constraint_dimensions){
        const const_dimension = constraint_dimensions[i];
        if(const_dimension.isNum){
            let constraint_object = { col: i, operator: "BETWEEN", val1: const_dimension.val[0], val2: const_dimension.val[1] };
            all_constraints.push(constraint_object);
        }else{
            let constraintarray = const_dimension.val;
            for(const j in constraintarray){
                let constraint_object = { col: i, operator: "EQUALS", val1: constraintarray[j], val2: constraintarray[j] };
                all_constraints.push(constraint_object);
            }
        
        }
    }
    
    ReactionCIMEBackendFromEnv.updatePOIConstraints(dataset.info.path, all_constraints).then((res) => {
        if(res.msg === "ok" && triggerDatasetUpdate != null){
            triggerDatasetUpdate({
                display: dataset.info.path,
                path: dataset.info.path,
                type: dataset.info.type,
                uploaded: true
            })
        }
    })
    
}