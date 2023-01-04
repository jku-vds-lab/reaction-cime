import { IVector } from 'projection-space-explorer';

export type ReactionVector = IVector & {
  circ_radius: number;
  hex: string;
};
