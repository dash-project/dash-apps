<?php

function getByName($name, $group) {
  foreach($group->content as $value) {
    if($name == $value->name) {
      return $value;
    }
  }
  return NULL;
}

$params_json = file_get_contents("php://input");
if(strlen($params_json) > 0) {
  $params = json_decode($params_json);
  $dimGroup = getByName("dim_group",$params);

  $pattern = " -s ".getByName("pattern",$params)->value->name;
  $size    = " -n ".getByName("size",$dimGroup)->value[0]." ".getByName("size",$dimGroup)->value[1];
  $units   = " -u ".getByName("units",$dimGroup)->value[0]." ".getByName("units",$dimGroup)->value[1];
  $tile    = " -t ".getByName("tile",$dimGroup)->value[0]." ".getByName("tile",$dimGroup)->value[1];
  $blocked = (getByName("blocked_display",$params)->value)?" -b":"";
  $balance = (getByName("balance_extents",$params)->value)?" -e":"";

  header("Content-Type: text/plain");
  //echo "./pattern-visualizer -p".$blocked.$pattern.$balance.$units.$size.$tile;
  passthru("./pattern-visualizer -p".$blocked.$pattern.$balance.$units.$size.$tile);
}
?>
