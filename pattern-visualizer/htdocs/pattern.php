<?php
$pattern = " -s ".$_GET['pattern'];
$size    = " -n ".$_GET['size_0']." ".$_GET['size_1'];
$units   = " -u ".$_GET['units_0']." ".$_GET['units_1'];
$tile    = " -t ".$_GET['tile_0']." ".$_GET['tile_1'];
$blocked = (strcmp($_GET['blocked_display'],"true") == 0)?" -b":"";
$balance = (strcmp($_GET['balance_extents'],"true") == 0)?" -e":"";

header("Content-Type: image/svg+xml");
passthru("./pattern-visualizer -p".$blocked.$pattern.$balance.$units.$size.$tile);
?>
