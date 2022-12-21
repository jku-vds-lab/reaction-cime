/**
 * Duck file for the LineUp input data
 */

import { filter } from 'lodash';

// const SET_DATA = "ducks/lineUpInput/SET_DATA"
// const SET_COLUMNS = "ducks/lineUpInput/SET_COLUMNS"
const SET_DUMP = 'ducks/lineUpInput/SET_DUMP';
const SET_FILTER = 'ducks/lineUpInput/SET_FILTER';
const UPDATE_FILTER = 'ducks/lineUpInput/UPDATE_FILTER';
const SET_LINEUP = 'ducks/lineUpInput/SET_LINEUP';
const SET_UPDATE = 'ducks/lineUpInput/SET_UPDATE';

// export const setLineUpInput_data = input => ({
//     type: SET_DATA,
//     input: input
// });

// export const setLineUpInput_columns = input => ({
//     type: SET_COLUMNS,
//     input: input
// });

export const setLineUpInputDump = (input) => ({
  type: SET_DUMP,
  input,
});

export const setLineUpInputFilter = (input) => ({
  type: SET_FILTER,
  input,
});

export const updateLineUpInputFilter = (input) => ({
  type: UPDATE_FILTER,
  input,
});

export const setLineUpInputLineup = (input) => ({
  type: SET_LINEUP,
  input,
});

export const setLineUpInputUpdate = (input) => ({
  type: SET_UPDATE,
  input,
});

const initialState: LineUpType = {
  // data: null,
  // columns: null,
  dump: '',
  filter: null,
  previousfilter: null,
  lineup: null,
  update: 0,
};
export type LineUpType = {
  // data: Vect[],
  // columns: [],
  dump: string;
  filter: object;
  previousfilter: object;
  lineup: any;
  update: number;
};

const lineUpInput = (state = initialState, action = undefined): LineUpType => {
  switch (action.type) {
    // case SET_DATA:
    //     return {...state, data: action.input}
    // case SET_COLUMNS:
    //     return {...state, columns: action.input}
    case SET_DUMP:
      return { ...state, dump: action.input };
    case SET_FILTER:
      return { ...state, previousfilter: { ...state.filter }, filter: action.input };
    case UPDATE_FILTER:
      if (state.filter && Object.keys(state.filter).includes(action.input.key)) {
        if (state.filter[action.input.key] === action.input.val_old) {
          const filterNew = { ...filter };
          filterNew[action.input.key] = action.input.val_new;
          return { ...state, filter: filterNew };
        }
      }
      return state;
    case SET_LINEUP:
      return { ...state, lineup: action.input };
    case SET_UPDATE:
      return { ...state, update: state.update + 1 };
    default:
      return state;
  }
};

export default lineUpInput;
