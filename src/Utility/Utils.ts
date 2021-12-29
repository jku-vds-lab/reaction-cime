
// https://stackoverflow.com/questions/31214677/download-a-reactjs-object-as-a-file
export const downloadImpl = (data: string, name: string, mimetype: string) => {
    var b = new Blob([data], { type: mimetype });
    var csvURL = window.URL.createObjectURL(b);
    let tempLink = document.createElement("a");
    tempLink.href = csvURL;
    tempLink.setAttribute("download", name);
    tempLink.click();
  };