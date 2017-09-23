// todo code cleaning
// todo callback function for more customisation
// todo concept graphics

var flexibleParams = new Object();

flexibleParams.config = new Object();
flexibleParams.config.maxQuantityDefault = 64;
flexibleParams.config.maxQuantitySoftlimit = 256;
flexibleParams.config.maxQuantityHardlimit = 4096;

flexibleParams.getIdSuffixDefault = function (quantityPosition) {
  var quantityIdSuffix = "";
  for(var i=0; i < quantityPosition.length; i++) {
    quantityIdSuffix += "_"+quantityPosition[i];
  }
  return quantityIdSuffix;
}
flexibleParams.getLabelSuffixDefault = function (quantityPosition) {
  var quantityLabelSuffix = "";
  for(var i=0; i < quantityPosition.length; i++) {
    if(i == 0) {
      quantityLabelSuffix += " "+quantityPosition[i];
    } else {
      quantityLabelSuffix += "."+quantityPosition[i];
    }
  }
  return quantityLabelSuffix;
}

// todo idea, introduce callback function to dynamically adjust those functions
/**
 * functions that allow to configure the numbering of labels
 * The default can be restored by setting them back to
 * flexibleParams.getIdSuffixDefault and flexibleParams.getLabelSuffixDefault.
 */
flexibleParams.config.getIdSuffix = flexibleParams.getIdSuffixDefault;
flexibleParams.config.getLabelSuffix = flexibleParams.getLabelSuffixDefault;

/**
 * stores the values recursive in the given argument param
 */
flexibleParams.storeParamValues = function (param,root,withIrrelevant,exceptPasswords,quantityCountList,quantityPosition) {
  var quantity = param.quantity;
  if(withIrrelevant == undefined) {
    withIrrelevant = false;
  }
  if(exceptPasswords == undefined) {
    exceptPasswords = false;
  }
  if(quantityCountList == undefined) {
    quantityCountList = [];
  }

  var quantityCount;
  if(quantity == undefined || (isNaN(quantity) && root.querySelector("#"+quantity) == null)) {
    quantityCount = 1;
  } else if(!isNaN(quantity)) {
    quantityCount = quantity;

    if(quantityCount > 1) {
      quantityCountList = quantityCountList.concat(quantityCount);
    }
  } else {
    quantityCount = flexibleParams.config.maxQuantityDefault;
    if(!isNaN(root.querySelector("#"+quantity).max)) {
      quantityCount = root.querySelector("#"+quantity).max*1;
    }
    if(!isNaN(param.maxQuantity) && param.maxQuantity*1 < quantityCount) {
      quantityCount = param.maxQuantity*1;
    }
    var tmpQuantityCount = quantityCount;

    if(!withIrrelevant && !isNaN(root.querySelector("#"+quantity).value)
                       && root.querySelector("#"+quantity).value > 0
                       && root.querySelector("#"+quantity).value < quantityCount
      ) {
      quantityCount = root.querySelector("#"+quantity).value*1;
    }

    if(tmpQuantityCount > 1) {
      quantityCountList = quantityCountList.concat(quantityCount);
    }
  }

  switch(param.type) {
    case "group":
      for(var i = 0; i < param.content.length; i++) {
        flexibleParams.storeParamValues(param.content[i],root,withIrrelevant,exceptPasswords,quantityCountList,quantityPosition);
      }
      break;
    case "password":
      if(exceptPasswords) {
        break;
      }
    case "text":
    case "hidden":
    case "number":
    case "range":
    case "logrange":
      param.value = flexibleParams.getParamValues(param,root,withIrrelevant,exceptPasswords,quantityCountList,quantityPosition);
      break;
    case "selection":
    case "selection-multiple":
    case "radio":
    case "checkbox":
    case "tab":
      param.value = flexibleParams.getParamValues(param,root,withIrrelevant,exceptPasswords,quantityCountList,quantityPosition);
      for(var i = 0; param.values != undefined && i < param.values.length; i++) {
        if(param.values[i].params != undefined &&
           (withIrrelevant || param.includeAll == true)
          ) {
          for(var j = 0; j < param.values[i].params.length; j++) {
            flexibleParams.storeParamValues(param.values[i].params[j],root,withIrrelevant,exceptPasswords,quantityCountList);
          }
        }
      }
      break;
    case "br":
    default:
      break;
  }
}

/**
 * returns an array with the values contained in the elements
 */
flexibleParams.getParamValues = function (param,root,withIrrelevant,exceptPasswords,quantityCountList,quantityPosition) {
  if(quantityPosition == undefined) {
    quantityPosition = [];
  }

  if(quantityCountList.length == 0) {
    var paramElem = root.querySelector("#"+param.name+flexibleParams.config.getIdSuffix(quantityPosition));
    if(param.type == "checkbox" && param.values == undefined) {
      return paramElem.checked;
    }
    if(param.type == "selection-multiple" || param.type == "selection" ||
       param.type == "radio" || param.type == "checkbox" || param.type == "tab") {
      var inpBased = param.type == "radio" || param.type == "checkbox" || param.type == "tab";
      var result = [];
      for(var i=0; i < param.values.length; i++) {
        var opt;
        if(!inpBased) {
          opt = paramElem.options[i];
        } else {
          opt = root.querySelector("#"+param.name+flexibleParams.config.getIdSuffix(quantityPosition)+param.values[i].name);
        }
        if((!inpBased && opt.selected) ||
           ( inpBased && opt.checked)) {
          var params = undefined;
          if(param.values[i].params != undefined && !withIrrelevant && !param.includeAll) {
            params = JSON.parse(JSON.stringify(param.values[i].params));
            var tmpGroup = {"type": "group", "content": params};
            flexibleParams.storeParamValues(tmpGroup,root,withIrrelevant,exceptPasswords,[],quantityPosition);
          }

          result.push({"name": opt.value, "params": params});
        }
      }
      if(param.type != "selection-multiple" && param.type != "checkbox") {
        return result[0];
      }
      return result;
    } else {
      return paramElem.value;
    }
  } else {
    var result = [];
    for(var i=0; i < quantityCountList[0]; i++) {
      result.push(flexibleParams.getParamValues(param,root,withIrrelevant,exceptPasswords,quantityCountList.slice(1),quantityPosition.concat(i)));
    }
    return result;
  }
}

/**
 * truncates the param object to reduce size
 * this creates a copy of the object as this process discards information
 */
flexibleParams.truncateParam = function (param,withIrrelevant) {
  var result = {};
  result.name = param.name;
  result.type = param.type;
  switch(param.type) {
    case "group":
      result.content = [];
      for(var i=0; i < param.content.length; i++) {
        var subParam = flexibleParams.truncateParam(param.content[i]);
        if(subParam.name != undefined) {
          result.content.push(subParam);
        }
      }
      break;
    case "text":
    case "password":
    case "hidden":
    case "number":
    case "range":
    case "logrange":
      result.value = JSON.parse(JSON.stringify(param.value));
      break;
    case "selection":
    case "selection-multiple":
    case "radio":
    case "checkbox":
    case "tab":
      result.value = flexibleParams.truncateSelectionValueParams(param.value);
      if(param.values != undefined && (withIrrelevant || param.includeAll)) {
        result.values = [];
        for(var i=0; i < param.values.length; i++) {
          var tmpValue = {"name": param.values[i].name};
          if(param.values[i].params != undefined) {
            var tmpGroup = {"type": "group", "content": param.values[i].params};
            tmpValue.params = flexibleParams.truncateParam(tmpGroup).content;
          }
          result.values.push(tmpValue);
        }
      }
      break;
    case "br":
      return {};
      break;
    default:
      console.warn("unknown parameter type: "+param.type);
      break;
  }
  return result;
}

flexibleParams.truncateSelectionValueParams = function (value) {
  if(Array.isArray(value)) {
    var result = [];
    for(var i=0; i < value.length; i++) {
      result.push(flexibleParams.truncateSelectionValueParams(value[i]));
    }
    return result;
  } else {
    if(value == true || value == false) {
      return value;
    }
    var result = {"name": value.name, "params": undefined};
    if(value.params != undefined) { // todo perhaps only return value.name otherwise
      var tmpGroup = {"type": "group", "content": value.params};
      result.params = flexibleParams.truncateParam(tmpGroup).content;
    }
    return result;
  }
}

/**
 * creates the specified param and calls the appropriate function
 * supported parameter types are:
 *  group
 *  text
 *  password
 *  number
 *  range
 *  logrange
 *  selection
 *  selection-multiple
 *  radio
 *  checkbox
 *  tab
 *  br - meta
 */
flexibleParams.createParam = function (param,quantityPosition,root,linebreak) {
  var quantity = param.quantity;
  var result = undefined;
  if(linebreak == undefined) {
    linebreak = {"val": false};
  }
  switch(param.type) {
    case "group":       // fieldset
      result = flexibleParams.createGroup(param,quantityPosition,root,linebreak);
      break;
    case "text":
    case "password":
    case "hidden":
    case "number":
    case "range":
    case "logrange":    // combined number + range
      if(quantity == undefined ||
         (isNaN(quantity) && root.querySelector("#"+quantity) == null)
        ) {
        result = flexibleParams.createInputParam(param,quantityPosition,linebreak);
        break;
      }/* else {
        // see else branch in next case
      }*/
    case "selection":   // implicit group for additional parameters
    case "selection-multiple":
    case "radio":
    case "checkbox":
    case "tab":
      if(quantity == undefined ||
         (isNaN(quantity) && root.querySelector("#"+quantity) == null)
        ) {
        result = flexibleParams.createSelection(param,quantityPosition,linebreak);
      } else {
        var tmpGroup = {"name":        param.name+"_quantityGroup",
                        "content":    [param],
                        "quantity":    param.quantity,
                        "maxQuantity": param.maxQuantity};
        param.quantity = undefined;
        result = flexibleParams.createGroup(tmpGroup,quantityPosition,root,linebreak);
        param.quantity = quantity;
      }
      break;
    case "br":    // css only linebreak
      linebreak.val = true;
      break;
    default:
      console.warn("unknown parameter type: "+param.type);
      break;
  }
  return result;
}

/**
 * creates a group (fieldset) with the specified contents and quantity
 */
flexibleParams.createGroup = function (param,quantityPosition,root,linebreak) {
  var quantity = param.quantity;
  if(quantityPosition == undefined) {
    quantityPosition = [];
  }
  if(linebreak == undefined) {
    linebreak = {"val": false};
  }

  var quantityIdSuffix = flexibleParams.config.getIdSuffix(quantityPosition);

  var group = document.createElement("fieldset");
  group.setAttribute("id",param.name+quantityIdSuffix);
  group.setAttribute("class","flexibleParams_container");
  if(linebreak.val) {
    group.setAttribute("class","flexibleParams_container"+" flexibleParams_linebreak");
    linebreak.val = false;
  }

  if(root == undefined) {
    root = group;
  }

  var quantityCount;
  var quantityDynamic = false;
  if(quantity == undefined || (isNaN(quantity) && root.querySelector("#"+quantity) == null)) {
    quantityCount = 1;
  } else if(!isNaN(quantity)) {
    quantityCount = quantity;
  } else {
    quantityDynamic = true;
    quantityCount = flexibleParams.config.maxQuantityDefault;
    if(!isNaN(root.querySelector("#"+quantity).max) && root.querySelector("#"+quantity).max != "") {
      quantityCount = root.querySelector("#"+quantity).max*1;
    }
    if(!isNaN(param.maxQuantity) && param.maxQuantity*1 < quantityCount) {
      quantityCount = param.maxQuantity*1;
    }

    var onchange = (function () {
      var tmpName = param.name;
      var tmpIdSuffix = quantityIdSuffix;
      return function (e) {
        if(!isNaN(e.target.value) && document.getElementById(tmpName+tmpIdSuffix+"_quantityHelper_"+(e.target.value-1)) != null) {
          document.getElementById(tmpName+tmpIdSuffix+"_quantityHelper_"+(e.target.value-1)).checked = true;
        };
      }
    })();

    root.querySelector("#"+quantity).addEventListener("change",onchange);
    root.querySelector("#"+quantity).addEventListener("input",onchange);
  }

  if(quantityCount > flexibleParams.config.maxQuantityHardlimit) {
    console.error("The given maximal quantity is too high: "+quantityCount+"\n"+
                  "You can adjust flexibleParams.config.maxQuantityHardlimit to change this behaviour.");
    quantityCount = flexibleParams.config.maxQuantitySoftlimit;
  } else if(quantityCount > flexibleParams.config.maxQuantitySoftlimit) {
    console.warn("The given maximal quantity is high: "+quantityCount+"\n"+
                 "You can supress this warning by changing flexibleParams.config.maxQuantitySoftlimit.");
  }

  for(var k = 0; k < quantityCount; k++) {
    if(k != 0 && param.content.length != 1) {
      linebreak.val = true;
    }
    var subgroupQuantityPosition = quantityPosition;
    if(quantityCount > 1) {
      subgroupQuantityPosition = subgroupQuantityPosition.concat(k);
    }

    // Main for loop that creates the content parameter
    for(var i = 0; i < param.content.length; i++) {
      var paramElem = flexibleParams.createParam(param.content[i],subgroupQuantityPosition,root,linebreak);
      if(paramElem != undefined) {
        group.appendChild(paramElem);
      }
    }

    if(quantityDynamic) {
      var opt_radio = document.createElement("input");
      opt_radio.setAttribute("type","radio");
      opt_radio.setAttribute("name",param.name+quantityIdSuffix+"_quantityHelper");
      opt_radio.setAttribute("id",param.name+quantityIdSuffix+"_quantityHelper_"+k);
      opt_radio.setAttribute("class","flexibleParams_quantityHelper");
      if(!isNaN(root.querySelector("#"+quantity).value) && (root.querySelector("#"+quantity).value*1)-1 == k) {
        opt_radio.checked = true;
      }

      group.appendChild(opt_radio);
    }
  }

  return group;
}

/**
 * creates a input field in a container with specified type and values
 */
flexibleParams.createInputParam = function (param, quantityPosition, linebreak) {
  var name = param.name;
  var type = param.type;
  var value = param.value;
  var min = param.min;
  var max = param.max;
  var step = param.step;
  if(quantityPosition == undefined) {
    quantityPosition = [];
  }

  var quantityIdSuffix = flexibleParams.config.getIdSuffix(quantityPosition);

  for(var i=0; i < quantityPosition.length; i++) {
    if(Array.isArray(value) && quantityPosition[i] < value.length) {
      value = value[quantityPosition[i]];
    } else {
      break;
    }
  }

  var cont = document.createElement("div");
  var input;
  var input_range;

  cont.setAttribute("class","flexibleParams_container");
  if(linebreak.val) {
    cont.setAttribute("class","flexibleParams_container"+" flexibleParams_linebreak");
    linebreak.val = false;
  }

  if(type == "logrange") {
    input = flexibleParams.createInput(name+quantityIdSuffix,"number",value,min,max,step);
    if(typeof logrange == "object") {
      input_range = flexibleParams.createInput(name+"_range"+quantityIdSuffix,"range",undefined,0,6,0.01);
      if(param.snappoints != undefined) {
        input_range.dataset.snappoints = JSON.stringify(param.snappoints);
      }

      input.addEventListener("input", function() { logrange.number_modified(input,input_range) });
      input_range.addEventListener("input", function() { logrange.range_modified(input,input_range) });
      input_range.addEventListener("change", function() { logrange.update_range(input,input_range,input_range.value) });
      logrange.number_modified(input,input_range);
    } else {
      console.error("include the script logrange.js before calling this script");
    }
  } else {
    input = flexibleParams.createInput(name+quantityIdSuffix,type,value,min,max,step);
  }
  if(param.label) {
    var label = document.createElement("label");
    var quantityLabelSuffix = flexibleParams.config.getLabelSuffix(quantityPosition);
    label.appendChild(document.createTextNode(param.label+quantityLabelSuffix));
    label.setAttribute("for",name+quantityIdSuffix);
    label.setAttribute("class","flexibleParams_label");

    cont.appendChild(label);
  }
  cont.appendChild(input);
  if(type == "logrange" && typeof logrange == "object") {
    cont.appendChild(input_range);
  }

  return cont;
}

flexibleParams.createInput = function (name,type,value,min,max,step) {
  var input = document.createElement("input");
  input.setAttribute("name",name);
  input.setAttribute("id",name);
  input.setAttribute("class","flexibleParams_input");
  input.setAttribute("type",type);
  if(value) {
    input.setAttribute("value",value);
  }
  if(min != undefined) {
    input.setAttribute("min",min);
  }
  if(max != undefined) {
    input.setAttribute("max",max);
  }
  if(step != undefined) {
    input.setAttribute("step",step);
  }

  return input;
}

/**
 * creates a selection field with specified values
 */
flexibleParams.createSelection = function (param, quantityPosition, linebreak) {
  var name = param.name;
  var type = param.type;
  var values = param.values;
  var size = param.size;
  var value = param.value;
  if(quantityPosition == undefined) {
    quantityPosition = [];
  }

  var multiple = false;
  if(type == "selection-multiple" || type == "checkbox") {
    multiple = true;
  }

  var quantityIdSuffix = flexibleParams.config.getIdSuffix(quantityPosition);

  for(var i=0; i < quantityPosition.length; i++) {
    if(Array.isArray(value) && quantityPosition[i] < value.length &&
       (!multiple ||
       ( multiple && Array.isArray(value[i])))
      ) {
      value = value[quantityPosition[i]];
    } else {
      break;
    }
  }

  var set = document.createElement("fieldset");
  var cont = document.createElement("div");

  set.setAttribute("class","flexibleParams_container");
  if(linebreak.val) {
    set.setAttribute("class","flexibleParams_container"+" flexibleParams_linebreak");
    linebreak.val = false;
  }
  cont.setAttribute("class","flexibleParams_container");

  if(param.label) {
    var label = document.createElement("label");
    var quantityLabelSuffix = flexibleParams.config.getLabelSuffix(quantityPosition);
    label.appendChild(document.createTextNode(param.label+quantityLabelSuffix));
    if(type == "selection" || type == "selection-multiple" || (type == "checkbox" && param.values == undefined)) {
      label.setAttribute("for",name+quantityIdSuffix);
    }
    label.setAttribute("class","flexibleParams_label");

    cont.appendChild(label);
  }
  if(type == "selection" || type == "selection-multiple") {
    var sel = document.createElement("select");
    sel.setAttribute("name",name+quantityIdSuffix);
    sel.setAttribute("id",name+quantityIdSuffix);
    sel.setAttribute("class","flexibleParams_input");
    if(type == "selection-multiple") {
      sel.multiple = true;
    }
    if(size != undefined) {
      sel.setAttribute("size",size);
    }
    cont.appendChild(sel);
  } else if(type == "checkbox" && param.values == undefined) {
    var opt = document.createElement("input");
    opt.setAttribute("name",name+quantityIdSuffix);
    opt.setAttribute("id",name+quantityIdSuffix);
    opt.setAttribute("type",type);
    opt.setAttribute("value","true");
    if(param.value == true) { // todo consider quantityPosition
      opt.checked = true;
    }
    cont.appendChild(opt);
  }

  set.appendChild(cont);

  for(var i=0; values != undefined && i < values.length; i++) {
    if(type == "selection" || type == "selection-multiple") {
      var opt = document.createElement("option");
      opt.setAttribute("value",values[i].name);
      if(values[i].label) {
        opt.appendChild(document.createTextNode(values[i].label));
      } else {
        opt.appendChild(document.createTextNode(values[i].name));
      }
      sel.appendChild(opt);
    } else if(type == "radio" || type == "checkbox" || type == "tab") {
      var opt = document.createElement("input");
      opt.setAttribute("name",name+quantityIdSuffix);
      if(type == "tab") {
        opt.setAttribute("type","radio");
        opt.setAttribute("class","flexibleParams_tabHelper");
        if(i == 0) {
          opt.checked = true;
        }
      } else {
        opt.setAttribute("type",type);
      }
      opt.setAttribute("value",values[i].name);
      opt.setAttribute("id",name+quantityIdSuffix+values[i].name);
      opt.addEventListener("change",function(e) {
        if(document.getElementById(e.target.id+"_selectionHelper") != null) {
          document.getElementById(e.target.id+"_selectionHelper").checked = e.target.checked;
        }
      });
      var optLabel = document.createElement("label");
      if(values[i].label) {
        optLabel.appendChild(document.createTextNode(values[i].label));
      } else {
        optLabel.appendChild(document.createTextNode(values[i].name));
      }
      optLabel.setAttribute("for",name+quantityIdSuffix+values[i].name);
      if(type == "tab") {
        optLabel.setAttribute("class","flexibleParams_tab")
      }
      cont.appendChild(opt);
      cont.appendChild(optLabel);
    }

    if(values[i].params != undefined) {
      var opt_radio = document.createElement("input");
      if(!multiple) {
        opt_radio.setAttribute("type","radio");
      } else {
        opt_radio.setAttribute("type","checkbox");
      }
      opt_radio.setAttribute("name",name+quantityIdSuffix+"_selectionHelper");
      opt_radio.setAttribute("id",name+quantityIdSuffix+values[i].name+"_selectionHelper");
      opt_radio.setAttribute("class","flexibleParams_selectionHelper");
      if(i == 0) {
        opt_radio.checked = true;
      }

      set.appendChild(opt_radio);
      var tmpGroup = {"name": name+quantityIdSuffix+"_SelectionGroup_"+i,
                      "content": values[i].params};
      var tmpI;
      if(value != undefined && typeof value == "object" &&
         value.name == values[i].name && value.params != undefined
        ) {
        tmpGroup.content = value.params;
      } else if(value != undefined && multiple &&
                -1 != (tmpI = value.findIndex(function(e) { // todo too bulky
                           if(typeof e != "object" || e == null) {
                               return false;
                           }
                           return e.name == values[i].name;
                       })) && value[tmpI].params != undefined
               ) {
        tmpGroup.content = value[tmpI].params;
      }
      set.appendChild(flexibleParams.createGroup(tmpGroup, quantityPosition));
    }

    if(value != undefined &&
       ((typeof value == "string" && value == values[i].name) ||
        (typeof value == "number" && value == i) ||
        (typeof value == "object" && value.name == values[i].name) ||
        (multiple && (value.includes(i) ||
                      value.includes(values[i].name) ||
                      -1 != value.findIndex(function(e) { // todo too bulky
                          if(typeof e != "object" || e == null) {
                              return false;
                          }
                          return e.name == values[i].name;
                      })
                     )
        )
       )
      ) {
      if(type == "selection" || type == "selection-multiple") {
        opt.selected = true;
      } else {
        opt.checked = true;
      }
      if(values[i].params != undefined) {
        opt_radio.checked = true;
      }
    }
  }

  if(type == "selection" || type == "selection-multiple") {
    sel.addEventListener("change",function(e) {
      for(var i = 0; i < e.target.options.length; i++) {
        if(document.getElementById(e.target.id+e.target.options[i].value+"_selectionHelper") != null) {
          document.getElementById(e.target.id+e.target.options[i].value+"_selectionHelper").checked = e.target.options[i].selected;
        }
      }
    });
  }

  return set;
}
