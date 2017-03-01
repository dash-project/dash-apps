<?php
header("Content-Type: image/svg+xml");
passthru("./pattern-visualizer -p -s ".$_GET['pattern']." -n ".$_GET['size_0']." ".$_GET['size_1']." -u ".$_GET['units_0']." ".$_GET['units_1']." -t ".$_GET['tile_0']." ".$_GET['tile_1']);
?>
