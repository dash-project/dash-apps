// todo custom round ?, logrange for comma values
var logrange = new Object();

logrange.config = new Object();
logrange.config.base = 2;
logrange.config.snap_dist = 0.08;

logrange.log = function (x) {
    switch(logrange.config.base) {
        case Math.E:
            return Math.log(x);
        case 2:
            if(Math.log2 != undefined) {
                return Math.log2(x);
            }
        case 10:
            if(Math.log10 != undefined && logrange.config.base == 10) {
                return Math.log10(x);
            }
        default:
            return Math.log(x)/Math.log(logrange.config.base);
    }
}

logrange.range_modified = function (number,range) {
    var snappoints = [];
    var snappoints_range = [];
    if(range.dataset.snappoint != undefined) {
        snappoints = JSON.parse(range.dataset.snappoints);
        snappoints.sort(function(a,b) {return a-b;});
        snappoints_range = snappoints.map(logrange.log);
    }

    var change = new InputEvent("change");

    for(var i=0; i < snappoints_range.length; i++) {
        if(Math.abs(snappoints_range[i] - range.value) < logrange.config.snap_dist) {
            number.value = snappoints[i];
            number.dispatchEvent(change);
            return;
        }
        if(range.value < snappoints_range[i] + logrange.config.snap_dist) {
            break;
        }
    }
    number.dispatchEvent(change);
    number.value = Math.round(Math.pow(logrange.config.base,range.value));
}

logrange.update_range = function (number,range,target) {

    var hardMax = Infinity;
    var hardMin = 0;

    if(!isNaN(number.max) && number.max != "") {
        hardMax = logrange.log(number.max*1);
    }
    if(!isNaN(number.min) && number.min != "" && number.min > 0) {
        hardMin = logrange.log(number.min*1);
    }


    while(target >= range.max*1 - 1 && range.max < hardMax) {
        range.min = Math.floor(range.min);

        range.min++;
        range.max++;

        if(range.max > hardMax) {
            range.max = hardMax;
        }
    }
    while(target <= range.min*1 + 1 && range.min > hardMin) {
        range.max = Math.ceil(range.max);

        range.min--;
        range.max--;

        if(range.min < hardMin) {
            range.min = hardMin;
        }
    }
}

logrange.number_modified = function (number,range) {
    var range_target = logrange.log(number.value*1);
    logrange.update_range(number,range,range_target);
    range.value = range_target;
}

