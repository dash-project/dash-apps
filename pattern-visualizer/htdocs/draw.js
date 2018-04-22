var pattern_combined;

function prepare_pattern(pattern) {
  pattern_combined = pattern;

  update_dims_list(pattern.dims,pattern.extents);

  // clear previous pattern

  // perhaps draw key
  update_key(pattern.maxUnits);

  // perhaps draw default view of pattern
  draw_pattern(pattern);
}

function update_dims_list(dims,extents) {
  var dims_list = document.getElementById("dims_list");
  while(dims_list.firstChild) {
    dims_list.removeChild(dims_list.firstChild);
  }

  for(var i=0; i < dims.length; i++) {
    var cont = document.createElement("div");
    var label = document.createElement("label");
    label.appendChild(document.createTextNode("Dimension "+dims[i]));
    var dimx = document.createElement("input");
    dimx.setAttribute("name","dimx_select");
    dimx.setAttribute("type","radio");
    dimx.setAttribute("class","dim_select");
    dimx.setAttribute("id","dimx_"+dims[i]);
    dimx.setAttribute("value",i);
    if(i == 1) {
      dimx.checked = true;
      document.getElementById("dimx").value = dims[i];
    }
    dimx.addEventListener("change",function(e) {
      var idx_dimx = e.target.value;
      var idx_dimy = document.getElementById("options").elements.namedItem("dimy_select").value;
      if(idx_dimx == idx_dimy) {
        var old_dimx = document.getElementById("dimx").value;
        var new_dimy = document.getElementById("dimy_"+old_dimx);
        if(new_dimy != undefined){
          new_dimy.checked = true;
          document.getElementById("dimy").value = old_dimx;
        }
      }
      document.getElementById("dimx").value = dims[e.target.value];
    });
    var dimx_label = document.createElement("label");
    dimx_label.setAttribute("for","dimx_"+dims[i]);
    dimx_label.appendChild(document.createTextNode("x"));
    var dimy = document.createElement("input");
    dimy.setAttribute("type","radio");
    dimy.setAttribute("class","dim_select");
    dimy.setAttribute("name","dimy_select");
    dimy.setAttribute("id","dimy_"+dims[i]);
    dimy.setAttribute("value",i);
    if(i == 0) {
      dimy.checked = true;
      document.getElementById("dimy").value = dims[i];
    }
    dimy.addEventListener("change",function(e) {
      var idx_dimx = document.getElementById("options").elements.namedItem("dimx_select").value;
      var idx_dimy = e.target.value;
      if(idx_dimx == idx_dimy) {
        var old_dimy = document.getElementById("dimy").value;
        var new_dimx = document.getElementById("dimx_"+old_dimy);
        if(new_dimx != undefined){
          new_dimx.checked = true;
          document.getElementById("dimx").value = old_dimy;
        }
      }
      document.getElementById("dimy").value = dims[e.target.value];
    });
    var dimy_label = document.createElement("label");
    dimy_label.setAttribute("for","dimy_"+dims[i]);
    dimy_label.appendChild(document.createTextNode("y"));
    var coord = document.createElement("input");
    coord.setAttribute("type","range");
    coord.setAttribute("id","coord_idx_"+i);
    coord.setAttribute("min",0);
    coord.setAttribute("max",extents[i]);
    coord.setAttribute("value",0);
    //coord.addEventListener("input",function(e){e.target.dataset.val = e.target.value});
    cont.appendChild(label);
    cont.appendChild(dimx);
    cont.appendChild(dimx_label);
    cont.appendChild(dimy);
    cont.appendChild(dimy_label);
    cont.appendChild(coord);
    dims_list.appendChild(cont);
  }
}

document.getElementById("dimx").addEventListener("change",function(e) {
  if(document.getElementById("dimx_"+e.target.value) != undefined) {
    var dimx = e.target.value;
    var dimy = document.getElementById("dimy").value;
    if(dimx != "" && dimx == dimy) {
      var old_idx_dimx = document.getElementById("options").elements.namedItem("dimx_select").value;
      var new_dimy = pattern_received.dims[old_idx_dimx];
      document.getElementById("dimy").value = new_dimy;
      document.getElementById("dimy_"+new_dimy).checked = true;
    }
    document.getElementById("dimx_"+e.target.value).checked = true;
  } else {
    e.target.value = "";
    if(pattern_received != undefined) {
      var idx_dimx = document.getElementById("options").elements.namedItem("dimx_select").value*1;
      e.target.value = pattern_received.dims[idx_dimx];
    }
  }
});
document.getElementById("dimy").addEventListener("change",function(e) {
  if(document.getElementById("dimy_"+e.target.value) != undefined) {
    var dimx = document.getElementById("dimx").value;
    var dimy = e.target.value;
    if(dimx != "" && dimx == dimy) {
      var old_idx_dimy = document.getElementById("options").elements.namedItem("dimy_select").value;
      var new_dimx = pattern_received.dims[old_idx_dimy];
      document.getElementById("dimx").value = new_dimx;
      document.getElementById("dimx_"+new_dimx).checked = true;
    }
    document.getElementById("dimy_"+e.target.value).checked = true;
  } else {
    e.target.value = "";
    if(pattern_received != undefined) {
      var idx_dimy = document.getElementById("options").elements.namedItem("dimy_select").value*1;
      e.target.value = pattern_received.dims[idx_dimy];
    }
  }
});

document.getElementById("blocked").addEventListener("change",function(e) {
  var nonBlockedStyle = document.getElementById("nonBlocked");
  if(nonBlockedStyle != undefined) {
    nonBlockedStyle.sheet.disabled = e.target.checked;
  }
});

function update_key(maxUnits) {
  var key = document.getElementById("key");
  while(key.firstChild) {
    key.removeChild(key.firstChild);
  }

  for(var i=0; i < maxUnits; i++) {
    var key_entry = document.createElement("div");
    var key_identifier = document.createElement("span");
    key_identifier.setAttribute("class","key_identifier u"+i);
    key_identifier.setAttribute("style","background-color: "+unit_color(i));
    var key_description = document.createTextNode("Unit "+i);
    key_entry.appendChild(key_identifier);
    key_entry.appendChild(key_description);
    key.appendChild(key_entry);
  }
}


/*var scale_value = 1.0;

function scale(e) {
  e.preventDefault();
  var prev_scale = scale_value;
  scale_value += 0.002*e.deltaY;
  if(prev_scale <= 1 && scale_value > 1) {
    document.getElementById("blocked").checked = true;
    draw_pattern(pattern_combined);
  } else if (prev_scale > 1 && scale_value <= 1) {
    document.getElementById("blocked").checked = false;
    draw_pattern(pattern_combined);
  }

  document.getElementById("result").style.transform = "scale("+scale_value+")";
}*/

function createElementSVG(elem) {
  return document.createElementNS("http://www.w3.org/2000/svg",elem);
}

var tile_size = 10;
var grid_base = tile_size + 2;
var fontsz    = 10;

function draw_pattern(pattern) {
  var coords = [];
  for(var i=0; i < pattern.dims.length; i++) {
    coords[i] = document.getElementById("coord_idx_"+i).value*1;
  }
  var idx_dimx = document.getElementById("options").elements.namedItem("dimx_select").value*1;
  var idx_dimy = document.getElementById("options").elements.namedItem("dimy_select").value*1;
  var blocked = document.getElementById("blocked").checked;

  var svg = createElementSVG("svg");

  // style, scripts
  create_defs(svg, pattern.maxUnits,blocked);

  // axes
  draw_axes(pattern, idx_dimx, idx_dimy, svg);

  var content = createElementSVG("g");
  content.setAttribute("transform","translate(24,24)");
  svg.appendChild(content);

  var svg_size = draw_blocks_regular(pattern, idx_dimx, idx_dimy, coords, content);

  if(pattern.memlayout != undefined) {
    draw_memlayout(pattern, idx_dimx, idx_dimy, coords, content);
  }

  var cont = document.getElementById("result");
  // clear previous result
  while(cont.firstChild) {
    cont.removeChild(cont.firstChild);
  }

  svg.setAttribute("width",svg_size["width"]+24+16);
  svg.setAttribute("height",svg_size["height"]+24+16);
  cont.appendChild(svg);

  if(blocked) {
    var nonBlockedStyle = document.getElementById("nonBlocked");
    if(nonBlockedStyle != undefined) {
      nonBlockedStyle.sheet.disabled = true;
    }
  }
}

function draw_axes(pattern, idx_dimx, idx_dimy, svg_elem) {
  var startx = 18;
  var starty = 18;
  var lenx = pattern.extents[idx_dimx] * grid_base + grid_base;
  var leny = pattern.extents[idx_dimy] * grid_base + grid_base;

  var axisX = createElementSVG("path");
  axisX.setAttribute("d","M"+startx+" "+starty+" h"+lenx);
  axisX.setAttribute("class","axis");
  axisX.setAttribute("marker-end","url(#arrowhead)");
  var axisX_label = createElementSVG("text");
  axisX_label.setAttribute("text-anchor","start");
  axisX_label.setAttribute("x",""+(startx + (lenx/3)));
  axisX_label.setAttribute("y",""+(starty - fontsz/2));
  axisX_label.setAttribute("font-size",""+fontsz);
  axisX_label.setAttribute("class","axis_label");
  axisX_label.appendChild(document.createTextNode("Dimension "+pattern.dims[idx_dimx]));
  var axisY = createElementSVG("path");
  axisY.setAttribute("d","M"+startx+" "+starty+" v"+leny);
  axisY.setAttribute("class","axis");
  axisY.setAttribute("marker-end","url(#arrowhead)");
  var axisY_label = createElementSVG("text");
  axisY_label.setAttribute("text-anchor","end");
  axisY_label.setAttribute("x",""+(startx - fontsz/2));
  axisY_label.setAttribute("y",""+(starty + (leny/3)));
  axisY_label.setAttribute("transform","rotate(-90,"+(startx-(fontsz/2))+","+(starty+(leny/3))+")");
  axisY_label.setAttribute("font-size",""+fontsz);
  axisY_label.setAttribute("class","axis_label");
  axisY_label.appendChild(document.createTextNode("Dimension "+pattern.dims[idx_dimy]));

  svg_elem.appendChild(axisX);
  svg_elem.appendChild(axisX_label);
  svg_elem.appendChild(axisY);
  svg_elem.appendChild(axisY_label);
}

function draw_blocks_regular(pattern, idx_dimx, idx_dimy, coords, svg_elem) {
  // strings for drawing in the right direction
  var dim_1, dim_2;
  var ext_1, ext_2;
  // position of dimension in depth of multi-dim-array
  var idx_1, idx_2;
  // position in current layer
  var pos_1, pos_2;
  // blocksize of current block
  var blocksize_1, blocksize_2;

  if(idx_dimx < idx_dimy) {
    idx_1 = idx_dimx;
    idx_2 = idx_dimy;
    dim_1 = "x";
    dim_2 = "y";
    ext_1 = "width";
    ext_2 = "height";
  } else if(idx_dimx > idx_dimy) {
    idx_1 = idx_dimy;
    idx_2 = idx_dimx;
    dim_1 = "y";
    dim_2 = "x";
    ext_1 = "height";
    ext_2 = "width";
  } else {
    console.warn("dimx == dimy");
    return;
  }

  pos_1 = 0;
  var slice_step_1 = get_slice_blocks_regular(0,idx_1,pattern.blocks,coords,pattern);
  for(var i = 0; i < slice_step_1.length; i++) {
    pos_2 = 0;
    var slice_step_2 = get_slice_blocks_regular(idx_1+1,idx_2,slice_step_1[i],coords,pattern);
    for(var j = 0; j < slice_step_2.length; j++) {
      var cur_block = get_slice_blocks_regular(idx_2+1,pattern.dims.length,slice_step_2[j],coords,pattern);
      blocksize_1 = pattern.blocksize[idx_1];
      blocksize_2 = pattern.blocksize[idx_2];
      if(cur_block.s) {
        blocksize_1 = cur_block.s[idx_1];
        blocksize_2 = cur_block.s[idx_2];
      }

      var group = createElementSVG("g");
      group.setAttribute("class","block_gr u"+cur_block.u);
      var block = createElementSVG("rect");
      block.setAttribute(dim_1, pos_1*grid_base);
      block.setAttribute(dim_2, pos_2*grid_base);
      block.setAttribute(ext_1, blocksize_1*grid_base - 2);
      block.setAttribute(ext_2, blocksize_2*grid_base - 2);
      block.setAttribute("class", "block");
      group.appendChild(block);

      var tiles = createElementSVG("g");
      tiles.setAttribute("class", "tile");
      for(var k=0; k < blocksize_2; k++) {
        for(var l=0; l < blocksize_1; l++) {
          var tile = createElementSVG("rect");
          tile.setAttribute(dim_1, pos_1*grid_base);
          tile.setAttribute(dim_2, pos_2*grid_base);
          tile.setAttribute("width",  tile_size);
          tile.setAttribute("height", tile_size);
          pos_1++;
          tiles.appendChild(tile);
        }
        pos_1 -= blocksize_1;
        pos_2++;
      }
      group.appendChild(tiles);
      svg_elem.appendChild(group);
    }
    pos_1 += blocksize_1;
  }

  var svg_size = new Array();
  svg_size[ext_1] = pos_1*grid_base;
  svg_size[ext_2] = pos_2*grid_base;
  return svg_size;
}

function get_first_block(blocks) {
  while(Array.isArray(blocks)) {
    blocks = blocks[0];
  }
  return blocks;
}

function get_slice_blocks_regular(idx_start, idx_dim, blocks, coords, pattern) {
  for(var i=idx_start; i < idx_dim; i++) {
    if(Array.isArray(blocks)) {
      var s=0;
      var target_coord = coords[i];
      var j=0;
      while(j < blocks.length) {
        var blocksize = pattern.blocksize[i];
        var first_block = get_first_block(blocks[j]);
        if(first_block.s) {
          blocksize = first_block.s[i];
        }
        s += blocksize;
        if(s > target_coord) {
          break;
        }
        j++;
      }
      if(j == blocks.length) {
        j = 0;
      }
      blocks = blocks[j];
    }
  }

  return blocks;
}

function draw_memlayout(pattern, idx_dimx, idx_dimy, coords, svg_elem) {
  for(var unit=0; unit < pattern.memlayout.length; unit++) {
    var local_memlayout = pattern.memlayout[unit];
    var mem_group = createElementSVG("g");
    mem_group.setAttribute("class","mem mem_u"+unit);

    var offset = 0;
    var contiguous = true;
    var startx = 0;
    var starty = 0;
    for(var i=0; i < local_memlayout.length; i++) {
      if(local_memlayout[i].o != undefined && offset != local_memlayout[i].o) {
        offset = local_memlayout[i].o;
        contiguous = false;
      }

      var pos = local_memlayout[i].p;
      if(match_slice(pos, coords, idx_dimx, idx_dimy)) {
        var endx = (pos[idx_dimx] * grid_base) + tile_size/2;
        var endy = (pos[idx_dimy] * grid_base) + tile_size/2;
        if(i != 0) {
          // draw line
          var line = createElementSVG("line");
          line.setAttribute("x1",startx);
          line.setAttribute("y1",starty);
          line.setAttribute("x2",endx);
          line.setAttribute("y2",endy);
          if(!contiguous) {
            line.setAttribute("class","discon");
          }
          mem_group.appendChild(line);
        }
        // draw point
        var circle = createElementSVG("circle");
        circle.setAttribute("cx",endx);
        circle.setAttribute("cy",endy);
        circle.setAttribute("r",1.5);
        mem_group.appendChild(circle);

        startx = endx;
        starty = endy;

        contiguous = true;
      } else {
        contiguous = false;
      }

      offset++;
    }
    svg_elem.appendChild(mem_group);
  }
}

function match_slice(coords1, coords2, idx_1, idx_2) {
  if(Array.isArray(coords1) && Array.isArray(coords2)) {
    for(var i=0; i < coords1.length && i < coords2.length; i++) {
      if(coords1[i] != coords2[i] && i != idx_1 && i != idx_2) {
        return false;
      }
    }
    return true;
  }
  return false;
}

// function createStyleSVG() {
  // var style = createElementSVG("style");
  // style.appendChild(document.createTextNode("/* "));
  // style.appendChild(document.createCDATASection(""));
  // style.appendChild(document.createTextNode(" */"));
  // return style;
// }*/

function create_defs(svg,maxUnits,blocked) {
  var defs = createElementSVG("defs");
  var unitStyle = createElementSVG("style");
  var unitStyle_content = "/* <![CDATA[ */\n";
  for(var i=0; i < maxUnits; i++) {
    unitStyle_content += ".u"+i+" { fill: "+unit_color(i)+" }\n";
    unitStyle_content += ".u"+i+":hover ~ .mem_u"+i+" { display: block }\n";
  }
  unitStyle_content += "/* ]]> */";
  unitStyle.appendChild(document.createTextNode(unitStyle_content));

  var defaultStyle = createElementSVG("style");
  var defaultStyle_content = "/* <![CDATA[ */\n";
  defaultStyle_content += ".axis { fill:grey;stroke:grey;stroke-width: 1 }\n";
  defaultStyle_content += ".axis_label { fill:grey;stroke-width: 0 }\n";
  defaultStyle_content += ".block { stroke-width: 0 }\n";
  defaultStyle_content += ".tile { stroke-width: 0 }\n";
  defaultStyle_content += ".mem { stroke: #E0E0E0; stroke-width: 1; fill: #E0E0E0; display: none }\n";
  defaultStyle_content += ".mem:hover { display: block }\n";
  defaultStyle_content += ".discon { stroke-dasharray: 4 }\n";
  defaultStyle_content += "/* ]]> */";
  defaultStyle.appendChild(document.createTextNode(defaultStyle_content));
  var blockedStyle = createElementSVG("style");
  var blockedStyle_content = "/* <![CDATA[ */\n";
  blockedStyle_content += ".block { display: block }\n";
  blockedStyle_content += ".tile { display: none }\n";
  blockedStyle_content += ".block_gr:hover > .block { opacity: 0 }\n";
  blockedStyle_content += ".block_gr:hover > .tile { display: block }\n";
  blockedStyle_content += "/* ]]> */";
  blockedStyle.appendChild(document.createTextNode(blockedStyle_content));
  var nonBlockedStyle = createElementSVG("style");
  nonBlockedStyle.setAttribute("id","nonBlocked");
  var nonBlockedStyle_content = "/* <![CDATA[ */\n";
  nonBlockedStyle_content += ".block { opacity: 0 }\n";
  nonBlockedStyle_content += ".tile { display: block }\n";
  nonBlockedStyle_content += "/* ]]> */";
  nonBlockedStyle.appendChild(document.createTextNode(nonBlockedStyle_content));

  var axisMarker = createElementSVG("marker");
  axisMarker.setAttribute("id","arrowhead");
  axisMarker.setAttribute("orient","auto");
  axisMarker.setAttribute("markerWidth","6");
  axisMarker.setAttribute("markerHeight","6");
  axisMarker.setAttribute("refX","0");
  axisMarker.setAttribute("refY","0");
  axisMarker.setAttribute("viewBox","-10 -15 30 30");
  var arrowHead = createElementSVG("path");
  arrowHead.setAttribute("d","M-10 -15 L20 0 L-10 15 L0 0 Z");
  arrowHead.setAttribute("class","axis");
  axisMarker.appendChild(arrowHead);

  defs.appendChild(unitStyle);
  defs.appendChild(defaultStyle);
  defs.appendChild(blockedStyle);
  defs.appendChild(nonBlockedStyle);
  defs.appendChild(axisMarker);
  svg.appendChild(defs);
}

function unit_color(unit) {
  var r = 0, g = 0, b = 0;
  switch (unit % 8) {
  case 0:
    r = 0x00;
    g = 0x72;
    b = 0xBD;
    break;
  case 1:
    r = 0xD9;
    g = 0x53;
    b = 0x19;
    break;
  case 2:
    r = 0xEB;
    g = 0xB1;
    b = 0x20;
    break;
  case 3:
    r = 0x7E;
    g = 0x2F;
    b = 0x8E;
    break;
  case 4:
    r = 0x77;
    g = 0xAC;
    b = 0x30;
    break;
  case 5:
    r = 0x4D;
    g = 0xBE;
    b = 0xEE;
    break;
  case 6:
    r = 0xA2;
    g = 0x14;
    b = 0x2F;
    break;
  case 7:
    r = 0x33;
    g = 0x6F;
    b = 0x45;
    break;
  }

  r += 20 * Math.floor(unit / 8);
  g += 20 * Math.floor(unit / 8);
  b += 20 * Math.floor(unit / 8);

  return "rgb("+(r%255)+","+(g%255)+","+(b%255)+")";
}
