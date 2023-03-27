import * as React from 'react';
import { Box, Button, Dialog, DialogActions, DialogContent, DialogTitle, IconButton } from '@mui/material';
import HelpIcon from '@mui/icons-material/Help';
import ReactMarkdown from 'react-markdown';
import readme from '../README.md?raw';

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
        <DialogTitle>Optional sizes</DialogTitle>
        <DialogContent>
          <ReactMarkdown>{readme}</ReactMarkdown>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleClose}>Close</Button>
        </DialogActions>
      </Dialog>
      <Box sx={{ position: 'absolute', right: 8, top: 8 }}>
        <IconButton size="small" onClick={handleClick}>
          <HelpIcon fontSize="inherit" />
        </IconButton>
      </Box>
    </>
  );
}
