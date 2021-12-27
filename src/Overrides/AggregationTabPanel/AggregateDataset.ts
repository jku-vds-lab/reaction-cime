import { FeatureType } from "projection-space-explorer";

export const DefaultFeatureLabel = "Default"

/**
 * Vector type.
 */
 export interface IVector {
    x: number;
    y: number;

}

type ColumnType = {
    distinct?: any
    metaInformation: any
    featureType: FeatureType
    range?: any
}

export class AggregateDataset {
    vectors: IVector[];
    bounds: { x; y; };
    columns: { [name: string] : ColumnType }

    constructor(vectors) {
        this.vectors = vectors;
        this.columns = {};

        this.bounds = {
            x: this.calculateRange("x"),
            y: this.calculateRange("y")
        };
        this.inferColumns();
    }

    getFeatureType(x) {
        if (typeof x === "number" || !isNaN(Number(x))) {
            return 'number'
        } else if ("" + new Date(x) !== "Invalid Date") {
            return 'date'
        } else {
            return 'arbitrary'
        }
    }

    inferColumns(){
        var header = Object.keys(this.vectors[0]);

        // Check for JSON header inside column, store it as key/value pair
        var metaInformation = header.reduce((map, value) => {
            let json = value.match(/[{].*[}]/)
            if (json != null) {
                let cutHeader = value.substring(0, value.length - json[0].length)

                this.vectors.forEach(vector => {
                    vector[cutHeader] = vector[value]
                    delete vector[value]
                })
                map[cutHeader] = JSON.parse(json[0])
            } else {
                map[value] = {"featureLabel": DefaultFeatureLabel}
            }
            return map
        }, {});

        // If data contains no x and y attributes, its invalid
        if (header.includes("x") && header.includes("y")) {
            this.vectors.forEach(vector => {
                vector.x = +vector.x;
                vector.y = +vector.y;
            });
        }else{
            console.log("You have to specify x and y columns in your dataset")//TODO: x and y are only needed for calculation of bounds -> set bounds to some default (e.g. [-1;1]) if x and y not given
        }
        //TODO: check, if the format of the dataset is ok (e.g. it has to be a cube, and only quantitative data is allowed...)

        Object.keys(metaInformation).forEach((key) => {
            // @ts-ignore
            this.columns[key] = { };

            this.retrieveFeatureTypes(metaInformation, key);

            //add some more information for each column depending on the feature type
            switch(this.columns[key].featureType){
                case FeatureType.Binary:
                    break;
                case FeatureType.Categorical:
                    this.columns[key].distinct = Array.from(new Set([...this.vectors.map(vector => vector[key])]));
                    break;
                case FeatureType.Date:
                    break;
                case FeatureType.Ordinal:
                    break;
                case FeatureType.Quantitative:
                    this.columns[key].range = this.calculateRange(key);
                    break;
                case FeatureType.String:
                    break;
                default:
                    break;
            }

        });

    }

    private retrieveFeatureTypes(metaInformation: {}, key: string) {
        const col_meta = metaInformation[key];
        if (col_meta?.dtype) {
            switch (col_meta.dtype) {
                case "numerical":
                    this.columns[key].featureType = FeatureType.Quantitative;
                    break;
                case "date":
                    this.columns[key].featureType = FeatureType.Date;
                    break;
                case "categorical":
                    this.columns[key].featureType = FeatureType.Categorical;
                    break;
                case "string":
                    this.columns[key].featureType = FeatureType.String;
                    break;
                default:
                    this.columns[key].featureType = FeatureType.String;
                    break;
            }
        } else if (metaInformation[key].range) {
            this.columns[key].range = metaInformation[key].range;
            this.columns[key].featureType = FeatureType.Quantitative;
        } else { //TODO: enhance this by automatic derivation of other featuretypes as well
            // infer for each feature whether it contains numeric, date, or arbitrary values
            var contains_number = {};
            var contains_date = {};
            var contains_arbitrary = {};
            this.vectors.forEach((r) => {
                const type = this.getFeatureType(r[key]);
                if (type === 'number') {
                    contains_number[key] = true;
                } else if (type === 'date') {
                    contains_date[key] = true;
                } else {
                    contains_arbitrary[key] = true;
                }
            });

            if (contains_number[key] && !contains_date[key] && !contains_arbitrary[key]) {
                // only numbers -> quantitative type
                this.columns[key].featureType = FeatureType.Quantitative;
            } else if (!contains_number[key] && contains_date[key] && !contains_arbitrary[key]) {
                // only date -> date type
                this.columns[key].featureType = FeatureType.Date;
            } else {
                // otherwise categorical
                this.columns[key].featureType = FeatureType.Categorical;
            }
        }
    }

    calculateRange(colName){
        var vector = this.vectors.map(vector => parseFloat(vector[colName])).filter(item => !isNaN(item));
        var min = Math.min(...vector);
        var max = Math.max(...vector);
        return {min: min, max: max};
    }

    map01(colName, vector){ // map the value to range [0;1]
        let div = this.columns[colName].range["max"]-this.columns[colName].range["min"];
        div = div > 0 ? div : 1;
        return (+vector[colName]-this.columns[colName].range["min"])/div;
    }

    normalize(colName, vector, lookup){ // normalize value to have 0 mean and unit sd; lookup is an object that may contain mean and std for certain columns
        let m, s;
                            
        if (colName in lookup) {
            m = lookup[colName].mean
            s = lookup[colName].std
        } else {
            m = mean(this.vectors.map(v => +v[colName]))
            s = std(this.vectors.map(v => +v[colName]))

            lookup[colName] = {
                mean: m,
                std: s
            }
        }

        if(s <= 0) // when all values are equal in a column, the standard deviation can be 0, which would lead to an error
            s = 1

        return {"result":(+vector[colName] - m) / s, "lookup":lookup};
    }

}

export function std(array) {
    const n = array.length
    const mean = array.reduce((a, b) => a + b) / n
    return Math.sqrt(array.map(x => Math.pow(x - mean, 2)).reduce((a, b) => a + b) / n)
}


export function mean(array) {
    const n = array.length
    const mean = array.reduce((a, b) => a + b) / n
    return mean
}