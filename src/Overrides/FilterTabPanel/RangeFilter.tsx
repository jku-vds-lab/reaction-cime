import { Box, Grid, IconButton, Slider, Tooltip, Typography } from "@mui/material";
import React from "react";
import DeleteIcon from '@mui/icons-material/Delete';
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { Dataset } from "projection-space-explorer";
import { formatLabelWithRange as formatLabel } from "../../Utility/Utils";


type Props = {
    col: string,
    value: number[],
    setValue: (value: number[]) => void,
    remove: (col:string) => void,
    dataset: Dataset,
};

export const RangeFilter = ({col, value, setValue, remove, dataset}:Props) => {
    const [min, setMin] = React.useState(value[0]);
    const [max, setMax] = React.useState(value[1]);

    React.useEffect(()=> {
        ReactionCIMEBackendFromEnv.loadValueRange(dataset.info.path, col).then((result) => { // have to load real possible range; TODO: could also include real range in the column info?
            setMin(result.min)
            setMax(result.max)
        })
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
                <Slider
                    getAriaLabel={() => col}
                    value={value}
                    onChange={handleChange}
                    valueLabelDisplay="auto"
                    aria-labelledby={"filter_"+col}
                    min={min}
                    max={max}
                    step={max == null || min == null || (Math.floor(max) === max || Math.floor(min) === min) ? 1 : (max-min)/100} // use step of 1 for natural numbers;
                    valueLabelFormat={(value) => formatLabel(value, min, max)}
                    marks={[{value: min ? min : 0, label: formatLabel(min ? min : 0, min, max)},{value: max ? max : 1, label: formatLabel(max ? max : 1, min, max)}]}
                />
            </Grid>
        </Grid>
    );
}