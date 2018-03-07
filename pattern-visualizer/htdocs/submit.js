function submitParams(e) {
  e.preventDefault();
  sendPOST();
  storeParams();
}

var pattern;

function sendPOST() {
  var req = new XMLHttpRequest();
  req.open("POST","pattern.php");
  req.addEventListener("load", function(e) {
    var responseText = e.target.responseText;
    console.log(responseText)
    var response = JSON.parse(responseText);
    if(response.success) {
      // update Options (dimensions / export link / save load share)
      pattern = response.pattern;
      draw_pattern(response.pattern);
      /*var svgText = responseText.slice(responseText.indexOf("}")+1,-1);
      var parser = new DOMParser();
      var result = parser.parseFromString(svgText,"image/svg+xml");
      var exportLink = document.createElement("a");
      exportLink.setAttribute("href","data:image/svg+xml;base64,"+window.btoa(svgText));
      exportLink.setAttribute("download",response.name);
      exportLink.setAttribute("target","_blank");
      exportLink.appendChild(document.createTextNode("Export as svg"));
      var imageDiv = document.createElement("div");
      imageDiv.appendChild(result.documentElement);
      displayResult([exportLink,imageDiv]);*/
    } else {
      displayResult([document.createTextNode(response.error)]);
    }
  });
  req.addEventListener("error", function(e) {
    var errorMsg = "Error loading visualized pattern, try again";
    console.warn(e);
    console.error(errorMsg);
    displayResult([document.createTextNode(errorMsg)]);
  });
  req.send(JSON.stringify(getParams()));
}

function displayResult(result) {
  var cont = document.getElementById("result");
  // clear previous result
  while(cont.firstChild) {
    cont.removeChild(cont.firstChild);
  }

  for(var i=0; i < result.length; i++) {
    cont.appendChild(result[i]);
  }
}

function createElementSVG(elem) {
  return document.createElementNS("http://www.w3.org/2000/svg",elem);
}

var scale_value = 1.0;

function scale(e) {
  e.preventDefault();
  var prev_scale = scale_value;
  scale_value += 0.002*e.deltaY;
  if(prev_scale <= 1 && scale_value > 1) {
    document.getElementById("blocked").checked = true;
    draw_pattern(pattern);
  } else if (prev_scale > 1 && scale_value <= 1) {
    document.getElementById("blocked").checked = false;
    draw_pattern(pattern);
  }

  document.getElementById("result").style.transform = "scale("+scale_value+")";
}

var tile_size = 10;
var grid_base = tile_size + 2;
var fontsz    = 10;

function draw_pattern(pattern) {
  var svg = createElementSVG("svg");
  svg.addEventListener("wheel",scale);
  var blocks = createElementSVG("g");
  var tiles = createElementSVG("g");
  if(document.getElementById("blocked").checked) {
    tiles.setAttribute("style","display: none");
  } else {
    blocks.setAttribute("style","display: none");
  }
  svg.appendChild(blocks);
  svg.appendChild(tiles);

  var posx, posy;
  var blocksize_x, blocksize_y;

  document.getElementById("dimx").innerHTML = "Dimension "+pattern.dims[0];
  document.getElementById("dimy").innerHTML = "Dimension "+pattern.dims[1];

  posx = 0;
  for(var i = 0; i < pattern.blocks.length; i++) {
    posy = 0;
    blocksize_x = pattern.blocksize[0];
    for(var j = 0; j < pattern.blocks[i].length; j++) {
      var cur_block = pattern.blocks[i][j];
      blocksize_x = pattern.blocksize[0];
      blocksize_y = pattern.blocksize[1];
      if(cur_block.s) {
        blocksize_x = cur_block.s[0];
        blocksize_y = cur_block.s[1];
      }

      var block = createElementSVG("rect");
      block.setAttribute("x", posx*grid_base);
      block.setAttribute("y", posy*grid_base);
      block.setAttribute("width",  blocksize_x*grid_base - 2);
      block.setAttribute("height", blocksize_y*grid_base - 2);
      block.setAttribute("style", "fill: "+unit_color(cur_block.u)+"; stroke-width: 0");
      blocks.appendChild(block);

      for(var k=0; k < blocksize_y; k++) {
        for(var l=0; l < blocksize_x; l++) {
          var tile = createElementSVG("rect");
          tile.setAttribute("x", posx*grid_base);
          tile.setAttribute("y", posy*grid_base);
          tile.setAttribute("width",  tile_size);
          tile.setAttribute("height", tile_size);
          tile.setAttribute("style", "fill: "+unit_color(cur_block.u)+"; stroke-width: 0");
          posx++;
          tiles.appendChild(tile);
        }
        posx -= blocksize_x;
        posy++;
      }
    }
    posx += blocksize_x;
  }

  svg.setAttribute("width",posx*grid_base);
  svg.setAttribute("height",posy*grid_base);

  var cont = document.getElementById("result");
  // clear previous result
  while(cont.firstChild) {
    cont.removeChild(cont.firstChild);
  }

  cont.appendChild(svg);
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
