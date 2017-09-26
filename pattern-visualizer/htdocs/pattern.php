<?php

function getByName($name, $group) {
  foreach($group->content as $value) {
    if($name == $value->name) {
      return $value;
    }
  }
  return NULL;
}

$params = json_decode(file_get_contents("php://input"));
$dimGroup = getByName("dim_group",$params);

$pattern = " -s ".getByName("pattern",$params)->value->name;
$size    = " -n ".getByName("size",$dimGroup)->value[0]." ".getByName("size",$dimGroup)->value[1];
$units   = " -u ".getByName("units",$dimGroup)->value[0]." ".getByName("units",$dimGroup)->value[1];
$tile    = " -t ".getByName("tile",$dimGroup)->value[0]." ".getByName("tile",$dimGroup)->value[1];
$blocked = (strcmp(getByName("blocked_display",$params)->value,"true") == 0)?" -b":"";
$balance = (strcmp(getByName("balance_extents",$params)->value,"true") == 0)?" -e":"";

header("Content-Type: image/svg+xml");
//echo "./pattern-visualizer -p".$blocked.$pattern.$balance.$units.$size.$tile;
passthru("./pattern-visualizer -p".$blocked.$pattern.$balance.$units.$size.$tile);
?>
