import { Grid } from "@mui/material";
import { DragAndDrop } from "projection-space-explorer";
import * as d3v5 from 'd3v5';

export var AggregatedDatasetDrop = ({ onChange }: { onChange(): void; }) => {
    return <Grid container item alignItems="stretch" justifyContent="center" direction="column" style={{ padding: '16px' }}>
        <DragAndDrop accept=".csv" handleDrop={(files) => {
            console.log("DragAndDrop")
            if (files == null || files.length <= 0) {
                return;
            }

            var file = files[0]
            var fileName = file.name as string

            var reader = new FileReader()
            reader.onload = (event) => {
                var content = event.target.result

                if (fileName.endsWith('csv')) {
                    console.log(content)
                    // console.log(d3v5.csvParse(content))
                    // new CSVLoader().resolveContent(content, onChange)
                }
            }

            reader.readAsText(file)
        }}>
            <div style={{ height: 200 }}></div>
        </DragAndDrop>
    </Grid>
}
