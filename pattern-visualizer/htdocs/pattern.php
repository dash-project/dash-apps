<?php
header("Content-Type: application/json");

$descriptorspec = array(
  0 => array("pipe", "r"),
  1 => array("pipe", "w")
);

$process = proc_open("./pattern-visualizer",$descriptorspec,$pipes);

if(is_resource($process)) {
  fwrite($pipes[0], file_get_contents("php://input"));
  fclose($pipes[0]);

  echo stream_get_contents($pipes[1]);
  fclose($pipes[1]);

  proc_close($process);
}
?>
