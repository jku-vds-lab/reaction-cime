import * as React from 'react';
import { Box, IconButton } from '@mui/material';
import HelpIcon from '@mui/icons-material/Help';

export function TabDocumentation({ value }: { value: number }) {
  const handleClick = () => {};

  return (
    <Box sx={{ position: 'absolute', right: 8, top: 8 }}>
      <IconButton size="small" onClick={handleClick}>
        <HelpIcon fontSize="inherit" />
      </IconButton>
    </Box>
  );
}
