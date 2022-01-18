import React from "react";
import ReactDOM from "react-dom";
import { ReactionCIMEApp } from "./ReactionCIMEApp";
import "./index.css";
import { version } from '../package.json';

/**
 * Main entry point for the app.
 */
ReactDOM.render(<ReactionCIMEApp />, document.getElementById("root"));
console.log("Reaction CIME version:", version)