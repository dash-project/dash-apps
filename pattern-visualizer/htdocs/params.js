var root;
var params;

function loadParams() {
  if(history.state != null) {
    params = history.state;
    appendParams();
  } else {
    var req = new XMLHttpRequest();
    req.open("GET","params.json");
    req.addEventListener("load", function(e) {
      paramsJSON = e.target.responseText;
      params = JSON.parse(paramsJSON);
      appendParams();
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

function appendParams() {
  root = flexibleParams.createParam(params);
  var end = document.getElementById("flexibleParams_end")
  document.getElementById("params").insertBefore(root,end);

  window.addEventListener("beforeunload",function(e) {
    if(params != undefined) {
      flexibleParams.storeParamValues(params,root,true,true);
      history.replaceState(params,"");
    }
  });

  window.addEventListener("popstate",function(e) {
    if(e.state != undefined) {
      params = e.state;
      var old_root = root;
      root = flexibleParams.createParam(params);
      document.getElementById("params").replaceChild(root,old_root);
    }
  });
}

function storeParams() {
  if(params != undefined) {
    flexibleParams.storeParamValues(params,root,true,true);
    history.replaceState(params,"");
    history.pushState(params,"");
  }
}

function getParams() {
  if(params != undefined) {
    flexibleParams.storeParamValues(params,root);
    return flexibleParams.truncateParam(params);
  }
}

