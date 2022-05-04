import { ReactionCIMEBackendFromEnv } from "../Backend/ReactionCIMEBackend";
import * as d3v5 from "d3v5";

export const PSE_BLUE = "#007dad" // "#1f77b4"
export const LIGHT_GREY = "#DDDDDD"


export const formatLabel = (value: number) => {
    if(value < 0.01){
        return d3v5.format(".2e")(value);
    }
    if(value >= 100){
        return Math.round(value)
    }
    
    return Math.round(value*100)/100;
}

export const formatLabelWithRange = (value, min, max) => {
    if(max == null || min == null)
        return Math.round(value)

    const step_size = (max-min)/100;
    if(step_size < 0.01){
        // https://github.com/d3/d3-format
        return d3v5.format(".2e")(value);
    }
    if(step_size >= 1)
        return Math.round(value)
    
    return Math.round(value*100)/100
}

export function arrayEquals(a, b) {
    return (
      Array.isArray(a) &&
      Array.isArray(b) &&
      a.length === b.length &&
      a.every((val, index) => val === b[index])
    );
  }

export function save_smiles_lookup_table(files:FileList){
    if (files == null || files.length <= 0) {
        return;
    }
    var file = files[0];

    const fileReader = new FileReader()
    fileReader.onload = (e) => {
        localStorage.setItem("smiles_lookup", e.target.result.toString());
        alert("Successfuly uploaded")
    }
    fileReader.readAsBinaryString(file)
}

export function map_smiles_to_shortname(smiles:string):string {
    const smiles_lookup_str = localStorage.getItem("smiles_lookup");
    if(smiles_lookup_str == null)
        return smiles
    const smiles_lookup = d3v5.csvParse(smiles_lookup_str) as Array<{smiles:string, shortname:string}>;
    return smiles_lookup.find((pair) => pair.smiles === smiles)?.shortname;
}

export function map_shortname_to_smiles(shortname:string):string {
    const smiles_lookup_str = localStorage.getItem("smiles_lookup");
    if(smiles_lookup_str == null)
        return shortname
    const smiles_lookup = d3v5.csvParse(smiles_lookup_str) as Array<{smiles:string, shortname:string}>;
    return smiles_lookup.find((pair) => pair.shortname === shortname)?.smiles;
}

/**
 * This is merely a helper function to decompose the code into smaller individual segments.
 * It is used to handle the changing background selection prop, which might, e.g., be triggered when a user clicks the k-nearest neighbor option in the context menu.
 * Specifically, it checks whether the parameters are correct, and if so, sends a query to the db to fetch the k-nearest entires to the click, with k being defined by a textfield.
 * The response triggers the download of a csv with these entries.
 * Afterwards, the background selection is reset, to make sure other prop updates do not trigger this db query and download.
 * @param {any} coords - The prop that holds x and y coordinates of the clicks
 * @returns {void} - no return value
 */
export function handleBackgroundSelectionDownload(coords: any, filename: string) {
    // if input for checking k-nearest neighbors (x,y coordinates and k) are not undefined
    if (
        typeof coords?.x !== "undefined" &&
        typeof coords?.y !== "undefined" &&
        (document.getElementById("knn-textfield") as HTMLInputElement)?.value !==
        "undefined"
    ) {
        let k = +(document.getElementById("knn-textfield") as HTMLInputElement)
        ?.value;
        // if input k is neither integer nor below 1
        if (k < 1 || k % 1 !== 0) {
            // warn user
            alert("Invalid input for k-nearest neighbors.");
        } else {
            // otherwise send request to db and download response in browser
            ReactionCIMEBackendFromEnv.getkNearestData(
                filename,
                coords?.x,
                coords?.y,
                (document.getElementById("knn-textfield") as HTMLInputElement)?.value
            );
        }
    }
}



/**
 * Creates a temporal dummy element and automatically triggers the download of the specified file with the specified content.
 * example usage: downloadImpl(JSON.stringify(response, null, 1), 'file.csv', 'text/csv')
 * // https://stackoverflow.com/questions/31214677/download-a-reactjs-object-as-a-file
 * @param {string} data - The content of the file to be downloaded. Make sure you call JSON.stringify(data) before passing it, if necessary.
 * @param {string} name - The name of the file including document file extension (such as "file.csv")
 * @param {string} mimetype - The mimetype of the file (such as "text/csv", see https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types)
 * @returns {void} - no return value
 */
export const downloadImpl = (data: string, name: string, mimetype: string) => {
    var b = new Blob([data], { type: mimetype });
    var csvURL = window.URL.createObjectURL(b);
    let tempLink = document.createElement("a");
    tempLink.href = csvURL;
    tempLink.setAttribute("download", name);
    tempLink.click();
    // TODO give created element a unique id and remove it again
  };



export const convert_to_rgb = (value: string | {r: number, g: number, b: number}):{r: number, g: number, b: number} => {
    if(Object.keys(value).includes("r") && Object.keys(value).includes("g") && Object.keys(value).includes("b"))
        return {"r": value["r"], "g": value["g"], "b": value["b"]};

    value = value.toString()
    if(value.startsWith("rgb")){
        var rgb = value.replace("rgb(", "")
        rgb = rgb.replace(")", "")
        rgb = rgb.replace(" ", "")
        var rgb_arr = rgb.split(",")
        return {"r": parseInt(rgb_arr[0]), "g": parseInt(rgb_arr[1]), "b": parseInt(rgb_arr[2])}
    }

    if(value.startsWith("#")){
        value = value.replace("#", "")
        var hex = value.match(/.{1,2}/g);
        return {"r": parseInt(hex[0], 16), "g": parseInt(hex[1], 16), "b": parseInt(hex[2], 16)}
    }
    
    console.log("error:", "format unknown ->", value)
    return {"r": 0, "g": 0, "b": 0};
}