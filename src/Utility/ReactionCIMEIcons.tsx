import * as React from 'react';
import SVG from 'react-inlinesvg';
import ReactionCIMEAggregateIcon from '../assets/pse-icon-aggregate.svg';
import ReactionCIMEFilterIcon from '../assets/pse-icon-filter.svg';

export const ReactionCIMEIcons = {
  Aggregate: () => <SVG src={ReactionCIMEAggregateIcon} />,
  Filter: () => <SVG src={ReactionCIMEFilterIcon} />, // TODO: change icon
};
