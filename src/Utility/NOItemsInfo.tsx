import { Typography } from '@mui/material';
import React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import { ReactionCIMEBackendFromEnv } from '../Backend/ReactionCIMEBackend';
import { AppState } from '../State/Store';

const mapStateToProps = (state: AppState) => ({
  dataset: state.dataset,
  globalLabels: state.globalLabels,
});

const mapDispatchToProps = (dispatch) => ({});

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>;

type Props = PropsFromRedux & {
  variant: 'all' | 'filterOutOfTotal';
};

export const NOItemsInfo = connector(({ dataset, globalLabels, variant }: Props) => {
  const [totalDataPoints, setTotalDataPoints] = React.useState(-1);
  const length = dataset?.vectors?.length ?? 0;
  React.useEffect(() => {
    if (dataset != null) {
      ReactionCIMEBackendFromEnv.loadNODatapoints(dataset.info.path).then((res) => {
        setTotalDataPoints(res.no_datapoints);
      });
    }
  }, [dataset]);

  let text;
  switch (variant) {
    case 'all':
      text = (
        <Typography color="textSecondary" variant="body2">
          Showing all <b>{totalDataPoints}</b> {globalLabels.itemLabelPlural} as aggregation
        </Typography>
      );
      break;
    case 'filterOutOfTotal':
      text = (
        <Typography color="textSecondary" variant="body2">
          Currently showing <b>{length}</b> out of <b>{totalDataPoints}</b> {globalLabels.itemLabelPlural} (~
          {Math.round((length / totalDataPoints) * 10000) / 100}% of the dataset).
        </Typography>
      );
      break;
    default:
      break;
  }

  return totalDataPoints >= 0 && text;
});
