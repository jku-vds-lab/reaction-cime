import { Table, TableBody, TableCell, TableHead, TableRow } from "@mui/material";
import { makeStyles } from "@mui/styles";
import { Dataset } from "projection-space-explorer";
import { connect, ConnectedProps } from "react-redux";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";
import { AppState } from "../../State/Store";


function genRows(aggregateSelection, aggregation, legendAttributes, dataset: Dataset) {
    // here is the code to retrieve the data for certain columns...
    ReactionCIMEBackendFromEnv.loadCategoryCount(dataset.info.path, "base_SMILES_index").then((data) => {
        console.log("count per category of base_SMILES_index for all datapoints")
        console.log(data)
    })

    ReactionCIMEBackendFromEnv.loadCategoryCountOfHex(dataset.info.path, "base_SMILES_index", aggregateSelection.x, aggregateSelection.y, aggregateSelection.circ_radius).then((data) => {
        console.log("count per category of base_SMILES_index for selected hexagon")
        console.log(data)
    })

    ReactionCIMEBackendFromEnv.loadDensity(dataset.info.path, "pred_mean_0").then((data) => {
        console.log("density data of pred_mean_0 for all datapoints")
        console.log(data)
    })

    ReactionCIMEBackendFromEnv.loadDensityOfHex(dataset.info.path, "pred_mean_0", aggregateSelection.x, aggregateSelection.y, aggregateSelection.circ_radius).then((data) => {
        console.log("density data of pred_mean_0 for selected hexagon")
        console.log(data)
    })
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
    // const rows = []
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