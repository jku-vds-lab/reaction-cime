import * as React from "react";
import { DatasetType, IVector, PSEPlugin } from "projection-space-explorer";
import { ReactionLegend } from "./ReactionLegend";

export class ReactionsPlugin extends PSEPlugin {
  type = DatasetType.Chem;

  createFingerprint(vectors: IVector[], scale: number, aggregate: boolean): JSX.Element {
      return <ReactionLegend selection={vectors} aggregate={aggregate} />;
  }
}