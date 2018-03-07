
function prepare_pattern(pattern) {
  update_dims_list(pattern.dims);

  // clear previous pattern

  // perhaps draw key

  // perhaps draw default view of pattern
  draw_pattern_blocks_regular(pattern);
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
    if(i == 0) {
      dimx.checked = true;
      document.getElementById("dimx").value = dims[i];
    }
    dimx.addEventListener("change",function(e) {
      var idx_dimx = e.target.value;
      var idx_dimy = document.getElementById("options").elements["dimy_select"].value;
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
    if(i == 1) {
      dimy.checked = true;
      document.getElementById("dimy").value = dims[i];
    }
    dimy.addEventListener("change",function(e) {
      var idx_dimx = document.getElementById("options").elements["dimx_select"].value;
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
    //coord.setAttribute("max",extents[i]);
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
      var old_idx_dimx = document.getElementById("options").elements["dimx_select"].value;
      var new_dimy = pattern_received.dims[old_idx_dimx];
      document.getElementById("dimy").value = new_dimy;
      document.getElementById("dimy_"+new_dimy).checked = true;
    }
    document.getElementById("dimx_"+e.target.value).checked = true;
  } else {
    e.target.value = "";
    if(pattern_received != undefined) {
      var idx_dimx = document.getElementById("options").elements["dimx_select"].value*1;
      e.target.value = pattern_received.dims[idx_dimx];
    }
  }
});
document.getElementById("dimy").addEventListener("change",function(e) {
  if(document.getElementById("dimy_"+e.target.value) != undefined) {
    var dimx = document.getElementById("dimx").value;
    var dimy = e.target.value;
    if(dimx != "" && dimx == dimy) {
      var old_idx_dimy = document.getElementById("options").elements["dimy_select"].value;
      var new_dimx = pattern_received.dims[old_idx_dimy];
      document.getElementById("dimx").value = new_dimx;
      document.getElementById("dimx_"+new_dimx).checked = true;
    }
    document.getElementById("dimy_"+e.target.value).checked = true;
  } else {
    e.target.value = "";
    if(pattern_received != undefined) {
      var idx_dimy = document.getElementById("options").elements["dimy_select"].value*1;
      e.target.value = pattern_received.dims[idx_dimy];
    }
  }
});

/*var scale_value = 1.0;

function scale(e) {
  e.preventDefault();
  var prev_scale = scale_value;
  scale_value += 0.002*e.deltaY;
  if(prev_scale <= 1 && scale_value > 1) {
    document.getElementById("blocked").checked = true;
    draw_pattern(pattern_received);
  } else if (prev_scale > 1 && scale_value <= 1) {
    document.getElementById("blocked").checked = false;
    draw_pattern(pattern_received);
  }

  document.getElementById("result").style.transform = "scale("+scale_value+")";
}*/

function createElementSVG(elem) {
  return document.createElementNS("http://www.w3.org/2000/svg",elem);
}

var tile_size = 10;
var grid_base = tile_size + 2;
var fontsz    = 10;

function draw_pattern_blocks_regular(pattern) {
  var svg = createElementSVG("svg");
  //svg.addEventListener("wheel",scale);
  var blocks = createElementSVG("g");
  var tiles = createElementSVG("g");
  if(document.getElementById("blocked").checked) {
    tiles.setAttribute("style","display: none");
  } else {
    blocks.setAttribute("style","display: none");
  }
  svg.appendChild(blocks);
  svg.appendChild(tiles);

  // strings for drawing in the right direction
  var dim_1, dim_2;
  var ext_1, ext_2;
  // position of dimension in depth of multi-dim-array
  var idx_1, idx_2;
  // position in current layer
  var pos_1, pos_2;
  // blocksize of current block
  var blocksize_1, blocksize_2;

  var idx_dimx = document.getElementById("options").elements["dimx_select"].value*1;
  var idx_dimy = document.getElementById("options").elements["dimy_select"].value*1;

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
  var slice_step_1 = get_slice(0,idx_1,pattern.blocks,pattern);
  for(var i = 0; i < slice_step_1.length; i++) {
    pos_2 = 0;
    var slice_step_2 = get_slice(idx_1+1,idx_2,slice_step_1[i],pattern);
    for(var j = 0; j < slice_step_2.length; j++) {
      var cur_block = get_slice(idx_2+1,pattern.dims.length,slice_step_2[j],pattern);
      blocksize_1 = pattern.blocksize[idx_1];
      blocksize_2 = pattern.blocksize[idx_2];
      if(cur_block.s) {
        blocksize_1 = cur_block.s[idx_1];
        blocksize_2 = cur_block.s[idx_2];
      }

      var block = createElementSVG("rect");
      block.setAttribute(dim_1, pos_1*grid_base);
      block.setAttribute(dim_2, pos_2*grid_base);
      block.setAttribute(ext_1, blocksize_1*grid_base - 2);
      block.setAttribute(ext_2, blocksize_2*grid_base - 2);
      block.setAttribute("style", "fill: "+unit_color(cur_block.u)+"; stroke-width: 0");
      blocks.appendChild(block);

      for(var k=0; k < blocksize_2; k++) {
        for(var l=0; l < blocksize_1; l++) {
          var tile = createElementSVG("rect");
          tile.setAttribute(dim_1, pos_1*grid_base);
          tile.setAttribute(dim_2, pos_2*grid_base);
          tile.setAttribute("width",  tile_size);
          tile.setAttribute("height", tile_size);
          tile.setAttribute("style", "fill: "+unit_color(cur_block.u)+"; stroke-width: 0");
          pos_1++;
          tiles.appendChild(tile);
        }
        pos_1 -= blocksize_1;
        pos_2++;
      }
    }
    pos_1 += blocksize_1;
  }

  svg.setAttribute(ext_1,pos_1*grid_base);
  svg.setAttribute(ext_2,pos_2*grid_base);

  var cont = document.getElementById("result");
  // clear previous result
  while(cont.firstChild) {
    cont.removeChild(cont.firstChild);
  }

  cont.appendChild(svg);
}

function get_first_block(blocks) {
  while(Array.isArray(blocks)) {
    blocks = blocks[0];
  }
  return blocks;
}

function get_slice_blocks_regular(idx_start, idx_dim, blocks, pattern) {
  for(var i=idx_start; i < idx_dim; i++) {
    if(Array.isArray(blocks)) {
      var s=0;
      var target_coord = document.getElementById("coord_idx_"+i).value*1;
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
