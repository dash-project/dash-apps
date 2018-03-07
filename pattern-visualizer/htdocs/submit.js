var params_sent;
var pattern_received;

function submitParams(e) {
  e.preventDefault();
  sendPOST();
  storeParams();
}

function sendPOST() {
  var req = new XMLHttpRequest();
  req.open("POST","pattern.php");
  req.addEventListener("load", function(e) {
    var responseText = e.target.responseText;
    console.log(responseText)
    var response = JSON.parse(responseText);
    if(response.success) {
      pattern_received = response.pattern;
      // update Options (dimensions / export link / save load share)
      prepare_pattern(response.pattern);
      /*var exportLink = document.createElement("a");
      exportLink.setAttribute("href","data:image/svg+xml;base64,"+window.btoa(svgText));
      exportLink.setAttribute("download",response.name);
      exportLink.setAttribute("target","_blank");
      exportLink.appendChild(document.createTextNode("Export as svg"));
      var imageDiv = document.createElement("div");
      imageDiv.appendChild(result.documentElement);
      displayResult([exportLink,imageDiv]);*/
    } else {
      var errorMsg = "An error occurred:\n\n"+response.error;
      console.warn(e);
      console.error(response.error);
      window.alert(errorMsg);
    }
  });
  req.addEventListener("error", function(e) {
    var errorMsg = "Error while loading the pattern from server, try again";
    console.warn(e);
    console.error(errorMsg);
    window.alert(errorMsg);
  });
  params_sent = JSON.stringify(getParams());
  req.send(params_sent);
}
