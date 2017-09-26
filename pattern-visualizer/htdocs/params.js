var root;
var params;

function loadParams() {
  if(history.state != null) {
    params = history.state;
    appendParams(params);
  } else {
    var req = new XMLHttpRequest();
    req.open("GET","params.json");
    req.addEventListener("load", function(e) {
      paramsJSON = e.target.responseText;
      params = JSON.parse(paramsJSON);
      appendParams(params);
    });
    req.addEventListener("error", function(e) {
      var errorMsg = "Error loading parameter, try reloading this page";
      console.warn(e);
      console.error(errorMsg);
      document.getElementById("params").appendChild(document.createTextNode(errorMsg));
    });
    req.send();
  }
}

function appendParams(params) {
  root = flexibleParams.createParam(params);
  document.getElementById("params").appendChild(root);

  window.addEventListener("beforeunload",function(e) {
    if(params != undefined) {
      flexibleParams.storeParamValues(params,root,true,true);
      history.replaceState(params,"");
    }
  });
}

function getParams() {
  flexibleParams.storeParamValues(params,root);
  return flexibleParams.truncateParam(params);
}

