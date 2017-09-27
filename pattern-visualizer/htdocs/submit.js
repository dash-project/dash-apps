function submitParams(e) {
  e.preventDefault();
  sendPOST();
}

function sendPOST() {
  var req = new XMLHttpRequest();
  req.open("POST","pattern.php");
  req.addEventListener("load", function(e) {
    var svgText = e.target.responseText;
    var parser = new DOMParser();
    var result = parser.parseFromString(svgText,"image/svg+xml");
    console.log(result);
    displayResult(result.documentElement);
  });
  req.addEventListener("error", function(e) {
    var errorMsg = "Error loading visualized pattern, try again";
    console.warn(e);
    console.error(errorMsg);
    displayResult(document.createTextNode(errorMsg));
  });
  req.send(JSON.stringify(getParams()));
}

function displayResult(result) {
  var cont = document.getElementById("result");
  // clear previous result
  while(cont.firstChild) {
    cont.removeChild(cont.firstChild);
  }

  cont.appendChild(result);
}

