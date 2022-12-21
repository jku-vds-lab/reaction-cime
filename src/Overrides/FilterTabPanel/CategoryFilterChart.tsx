import { Grid, IconButton, ToggleButton, Tooltip, Typography } from "@mui/material";
import React from "react";
import DeleteIcon from '@mui/icons-material/Delete';
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { Dataset } from "projection-space-explorer";
import { ReactBarChart } from "./D3Helpers/BarChart";


type Props = {
    col: string,
    value: string[],
    setValue: (value: string[]) => void,
    remove: (col:string) => void,
    dataset: Dataset,
};

export const CategoryFilterChart = ({col, value, setValue, remove, dataset}:Props) => {

    const [catValuesSample, setCatValuesSample] = React.useState(null)
    // const [catValuesTotal, setCatValuesTotal] = React.useState(null)

    React.useEffect(()=> {
        // load distinct features and their count in the currently filtered points -> better because it shows where users can still filter
        // problem: how to add filtered datavalues? -> maybe some outline of total values?
        ReactionCIMEBackendFromEnv.loadCategoryValues(dataset.info.path, col).then((result) => {
            const data = result.values.map((row) => {
                return {"feature": row, "count": dataset.vectors.filter((vector) => vector[col] === row).length}
            })
            setCatValuesSample(data)
        })
        // load distinct features and their total count in database -> does not make too much sense because the feature values are usually equally distributed
        // maybe as additional "context" information?
        // ReactionCIMEBackendFromEnv.loadCategoryCount(dataset.info.path, col).then((result) => {
            
        //     const data = result.map((row) => { 
        //         return {"feature": row[col], "count": row["count"]}
        //     })
        //     setCatValuesTotal(data)
        // })
    }, [col, dataset.info.path]);
    
    return (
        <Grid container paddingTop={0}>
            <Grid item xs={3} textAlign={"right"}>
                <Tooltip title={"Remove " + col + " filter"}>
                    <IconButton onClick={()=>{remove(col)}}><DeleteIcon fontSize="large"/></IconButton>
                </Tooltip>
            </Grid>
            <Grid item xs={9}>
                <Typography id={"filter_"+col} marginBottom={"0px"}>
                    {col}
                </Typography>
                <div>
                    {catValuesSample && <ReactBarChart data={catValuesSample} value={value} setValue={setValue} isSmiles={dataset.columns[col]?.metaInformation.imgSmiles}></ReactBarChart>}
                </div>
            </Grid>
        </Grid>
    );
}

