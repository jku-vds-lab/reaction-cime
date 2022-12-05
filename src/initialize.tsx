// import { App } from './app';

// eslint-disable-next-line @typescript-eslint/no-unused-vars
// const unused = new App({
//   name: 'CIME: ChemInformatics Model Explorer',
//   showAboutLink: false,
//   showHelpLink: false,
//   showOptionsLink: false,
//   showProvenanceMenu: false,
//   showReportBugLink: false,
//   showResearchDisclaimer: false,
//   enableProvenanceUrlTracking: false,
// });

import * as React from "react";
import ReactDOM from "react-dom";
import { ReactionCIMEApp } from "./ReactionCIMEApp";
import "./index.css";

ReactDOM.render(<ReactionCIMEApp />, document.getElementById("app"));