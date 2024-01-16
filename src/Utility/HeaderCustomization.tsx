import * as React from 'react';
import { Group, Anchor, Text, Stack } from '@mantine/core';
import bayerLogo from '../assets/bayer_logo_dark_text.svg';
import vdsLogo from '../assets/jku-vds-lab-logo.svg';

export function BuildInfoContent() {
  return (
    <Stack gap="xl">
      <Text>
        CIME4R is a web-based tool for interactive exploration of chemical reaction optimization data. CIME4R allows chemists to better understand the decision
        process of AI models and enhances human-AI collaboration.
      </Text>
    </Stack>
  );
}

export function BuildInfoLogos() {
  return (
    <Group align="center" wrap="nowrap">
      <Anchor href="https://bayer.com/" rel="noreferrer noopener" target="_blank" style={{ transform: 'scale(1.3)' }}>
        <img src={bayerLogo} alt="logo" style={{ height: '24px' }} />
      </Anchor>
      <Anchor href="https://jku-vds-lab.at/" rel="noreferrer noopener" target="_blank">
        <img src={vdsLogo} alt="logo" style={{ height: '24px' }} />
      </Anchor>
    </Group>
  );
}
