const SET = "ducks/aggregateColor/SET"

const aggregateColor = (state = null, action) => {
    switch (action.type) {
        case SET:
            return action.aggregateColor
        default:
            return state
    }
}


export const setAggregateColor = aggregateColor => ({
    type: SET,
    aggregateColor: aggregateColor
})

export default aggregateColor