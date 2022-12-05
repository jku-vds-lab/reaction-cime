import * as React from 'react';
import { DefaultLegend, IVector } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";
import { AggregateLegend } from "./AggregateLegend";
import { FeatureLegend } from "./FeatureLegend";

const mapStateToProps = (state: AppState) => ({
    aggregateSelection: state.selection?.currentAggregateSelection,
})
const mapDispatchToProps = (dispatch: any) => ({
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    aggregate: boolean;
    selection: IVector[];
}

export const ReactionLegend = connector(({ aggregate, selection, aggregateSelection }: Props) => {
    
    if(selection.length > 0){
        return <FeatureLegend selection={selection} aggregate={aggregate} />;
    }
    else if(aggregateSelection != null){
        return <AggregateLegend aggregate={aggregate} aggregateSelection={aggregateSelection}></AggregateLegend>
    }else{
        return <DefaultLegend></DefaultLegend>
    }
});