import { Table, TableBody, TableCell, TableHead, TableRow } from "@mui/material";
import { makeStyles } from "@mui/styles";
import { DefaultLegend, IVector } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { AppState } from "../../State/Store";


function genRows(vectors, aggregation, legendAttributes, dataset) {
    return [];
}

const mapStateToProps = (state: AppState) => ({
    // aggregateSelection: state.selection?.currentAggregateSelection,
    legendAttributes: state.genericFingerprintAttributes,
    dataset: state.dataset,
})
const mapDispatchToProps = (dispatch: any) => ({
})

const connector = connect(mapStateToProps, mapDispatchToProps);
type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
    aggregate: boolean;
    aggregateSelection: any;
}

const useStyles = makeStyles({
    table: {
      maxWidth: 288,
    },
    tableRow: {
      height: '66px',
    },
  });

export const AggregateLegend = connector(({ aggregate, aggregateSelection, legendAttributes, dataset }: Props) => {
    const classes = useStyles();
    const rows = genRows(aggregateSelection, aggregate, legendAttributes, dataset);

    return <div style={{ width: '100%', maxHeight: '100%', overflowY: 'scroll' }}>
        <div
        style={{
            width: '100%',
            // overflow: "auto"
        }}
        >
        <Table className={classes.table} aria-label="simple table" size="small">
            <TableHead />
            <TableBody>
            {rows.map((row) => (
                <TableRow className={classes.tableRow} key={row.feature}>
                    <TableCell component="th" scope="row">
                        <div style={{ maxWidth: 200 }}>
                        {row.feature}
                        <br />
                        <b>{row.category}</b>
                        </div>
                    </TableCell>
                    <TableCell>{row.char}</TableCell>
                </TableRow>
            ))}
            </TableBody>
        </Table>
        </div>
    </div>
});