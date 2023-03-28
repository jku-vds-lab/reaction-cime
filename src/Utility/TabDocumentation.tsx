import * as React from 'react';
import { Box, Button, Dialog, DialogActions, DialogContent, DialogTitle, IconButton, Tooltip, Typography } from '@mui/material';
import HelpIcon from '@mui/icons-material/Help';
import ReactMarkdown from 'react-markdown';
import DatasetReadme from '../Readme/Dataset_README.md?raw';
import AggregateReadme from '../Readme/Aggregate_README.md?raw';
import EncodingReadme from '../Readme/Encoding_README.md?raw';
import ProjectionReadme from '../Readme/Projection_README.md?raw';
import SelectionReadme from '../Readme/Selection_README.md?raw';
import TabularReadme from '../Readme/Tabular_README.md?raw';
import FilterReadme from '../Readme/Filter_README.md?raw';
import GroupReadme from '../Readme/Group_README.md?raw';

const readmeMap = {
  0: DatasetReadme,
  1: ProjectionReadme,
  2: EncodingReadme,
  3: GroupReadme,
  4: SelectionReadme,
  5: TabularReadme,
  6: FilterReadme,
  7: AggregateReadme,
};

export function TabDocumentation({ value }: { value: number }) {
  const [open, setOpen] = React.useState(false);

  const handleClick = () => {
    setOpen(true);
  };

  const handleClose = () => {
    setOpen(false);
  };

  return (
    <>
      <Dialog open={open} onClose={handleClose} maxWidth="xl">
        <DialogContent>
          <Typography>
            <ReactMarkdown>{readmeMap[value]}</ReactMarkdown>
          </Typography>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleClose}>Close</Button>
        </DialogActions>
      </Dialog>
      <Box sx={{ position: 'absolute', right: 8, top: 8 }}>
        <Tooltip title="Shows the documentation for this tab panel">
          <IconButton size="small" onClick={handleClick}>
            <HelpIcon fontSize="inherit" />
          </IconButton>
        </Tooltip>
      </Box>
    </>
  );
}
