import { Grid } from "@mui/material";
import * as d3v5 from 'd3v5';
import { DragAndDrop } from "projection-space-explorer";
import { AggregateDataset } from "./AggregateDataset";
// import DragAndDrop from "./AggregatedDSDragAndDrop";
import { ReactionCIMEBackendFromEnv } from "../../Backend/ReactionCIMEBackend";

export var AggregatedDatasetDrop = ({ onChange }: { onChange(dataset: AggregateDataset): void; }) => {
    return <Grid container item alignItems="stretch" justifyContent="center" direction="column" style={{ padding: '16px' }}>
        <DragAndDrop accept=".csv" handleDrop={(files) => {
            if (files == null || files.length <= 0) {
                return;
            }

            var file = files[0]
            var fileName = file.name as string

            var reader = new FileReader()
            reader.onload = (event) => {
                var content = event.target.result

                if (fileName.endsWith('csv')) {
                    onChange(d3v5.csvParse(content))
                }
            }

            reader.readAsText(file)
        }}>
            <div style={{ height: 200 }}></div>
        </DragAndDrop>
    </Grid>
}
