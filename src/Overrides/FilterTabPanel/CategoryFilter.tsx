import { Box, Grid, IconButton, Tooltip, Typography } from "@mui/material";
import React from "react";
import DeleteIcon from '@mui/icons-material/Delete';
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { Dataset } from "projection-space-explorer";
import * as d3v5 from "d3v5";


type Props = {
    col: string,
    value: number[],
    setValue: (value: number[]) => void,
    remove: (col:string) => void,
    dataset: Dataset,
};

export const CategoryFilter = ({col, value, setValue, remove, dataset}:Props) => {

    React.useEffect(()=> {

    }, [])
    
    const handleChange = (event, newValue) => {
        setValue(newValue);
    };
    return (
        <Grid container paddingTop={0}>
            <Grid item xs={3} textAlign={"right"}>
                <Tooltip title={"Remove " + col + " filter"}>
                    <IconButton onClick={()=>{remove(col)}}><DeleteIcon fontSize="large"/></IconButton>
                </Tooltip>
            </Grid>
            <Grid item xs={9}>
                <Typography id={"filter_"+col} marginBottom={"-5px"}>
                    {col}
                </Typography>
                <div>
                    TODO
                </div>
            </Grid>
        </Grid>
    );
}