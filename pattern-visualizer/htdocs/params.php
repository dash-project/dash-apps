<?php
header("Content-Type: application/json");
passthru("./pattern-visualizer -p");
?>
