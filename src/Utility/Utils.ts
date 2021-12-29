

/**
 * Creates a temporal dummy element and automatically triggers the download of the specified file with the specified content.
 * example usage: downloadImpl(JSON.stringify(response, null, 1), 'file.csv', 'text/csv')
 * // https://stackoverflow.com/questions/31214677/download-a-reactjs-object-as-a-file
 * @param {string} data - The content of the file to be downloaded. Make sure you call JSON.stringify(data) before passing it, if necessary.
 * @param {string} name - The name of the file including document file extension (such as "file.csv")
 * @param {string} mimetype - The mimetype of the file (such as "text/csv", see https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types)
 * @returns {void} no return value
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