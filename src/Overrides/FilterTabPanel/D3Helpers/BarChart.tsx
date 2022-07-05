import Tooltip, { tooltipClasses } from '@mui/material/Tooltip';
import * as d3 from "d3v5";
import React from "react";
import { ReactionCIMEBackendFromEnv } from "../../../Backend/ReactionCIMEBackend";
import { LIGHT_GREY, map_smiles_to_shortname, PSE_BLUE } from "../../../Utility/Utils";
import { styled } from '@mui/material/styles';

type Props = {
    data: {feature: string, count: number}[],
    value: string[],
    isSmiles: boolean,
    setValue: (value: string[]) => void,
};

export const ReactBarChart = ({data, value, isSmiles, setValue}: Props) => {
    const rel_width = 2;
    const rel_height = 1;

    const [xScale, setXScale] = React.useState<any>()
    const [yScale, setYScale] = React.useState<any>()

    React.useEffect(() => {
        if(data != null){
            console.log(data)
            // define x and y scales
            setXScale(() => d3.scaleBand()
                .domain(
                    data.map(function (d) { 
                        return d.feature; 
                    })
                )
                .range([0, rel_width]));
    
            setYScale(() => d3.scaleLinear()
                .domain([0,
                    d3.max(data, function (d) { 
                        return d.count; 
                    }), 
                ])
                .range([0, rel_height]));
        }

    }, [data])
    return <svg preserveAspectRatio="xMidYMid meet" width="100%" height="100" viewBox={`0 0 ${rel_width} ${rel_height}`}>
        {xScale && yScale &&
            <g>
                {
                    data.map((d) => {
                        return <LightTooltip followCursor={true} placement="right-start" style={{backgroundColor: "white"}} title={
                            <FeatureTooltipTitle feature={d.feature} isSmiles={isSmiles} count={d.count}></FeatureTooltipTitle>
                        }>
                            <rect 
                                x={xScale(d.feature)} 
                                y={rel_height-yScale(d.count)} 
                                width={xScale.bandwidth()} 
                                height={yScale(d.count)} 
                                style={{fill: value.includes(d.feature) ? PSE_BLUE : LIGHT_GREY, strokeWidth: xScale.bandwidth()*0.05, stroke: "white", cursor: "pointer"}}
                                onClick={() => {
                                    if(value.includes(d.feature)){ // if previously checked, remove from list
                                        setValue([...new Set(value.filter((v) => v !== d.feature))])
                                    }else{ // if previously unchecked, add to list
                                        setValue([...new Set([...value, d.feature])])
                                    }
                                }}
                            ></rect>
                        </LightTooltip>
                    })
                }
            </g>
        }
    </svg>
}

const LightTooltip = styled(({ className, ...props }: any) => (
    <Tooltip {...props} classes={{ popper: className }} />
  ))(({ theme }) => ({
    [`& .${tooltipClasses.tooltip}`]: {
      backgroundColor: theme.palette.common.white,
      color: 'rgba(0, 0, 0, 0.87)',
      boxShadow: theme.shadows[1],
      fontSize: 11,
    },
  }));

export const FeatureTooltipTitle = ({feature, count, isSmiles}) => {

    const [smilesImg, setSmilesImg] = React.useState<string>();

    React.useEffect(() => {
        if(isSmiles){
            let smiles = feature;
            ReactionCIMEBackendFromEnv.getStructureFromSmiles(smiles, false, null).then((x) => {
                if (x && x.length > 100) {
                    // check if it is actually long enogh to be an img
                    setSmilesImg(`url('data:image/jpg;base64,${x}')`);
                }else{
                    setSmilesImg(null)
                }
            });
        }
        
    }, [feature])

    return <React.Fragment>
        <div>Count: {count}</div>
        <div>Feature: {feature}</div>
        {isSmiles && <div>Short name: {map_smiles_to_shortname(feature)}</div>}
        {smilesImg && <div id={"smiles_${value.category}"} style={{width:"100%", height:"100px", backgroundSize: "contain", backgroundPosition: "center", backgroundRepeat: "no-repeat", backgroundImage: smilesImg}}></div>}
    </React.Fragment>
}