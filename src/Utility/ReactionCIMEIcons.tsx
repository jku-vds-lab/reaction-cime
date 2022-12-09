import ReactionCIMEAggregateIcon from '../assets/pse-icon-aggregate.svg'
import ReactionCIMEFilterIcon from '../assets/pse-icon-filter.svg'
import ReactionCIMETableIcon from '../assets/pse-icon-table.svg'
import SVG from 'react-inlinesvg';

export const ReactionCIMEIcons = {
    Table: () => <SVG src={ReactionCIMETableIcon}/>,
    Aggregate: () => <SVG src={ReactionCIMEAggregateIcon} />,
    Filter: () => <SVG src={ReactionCIMEFilterIcon}/>,
}