var ex01 = // todo complete
{
"name": "ex01",
"type": "group",
"content": [
{"name":"numUnits","label":"Number of Units","type":"number","min":1},
{"name":"memArrangement","label":"Memory Arrangement","type":"selection","values":[{"name":"ROW_MAJOR"},{"name":"COL_MAJOR"}]},
 {"type":"br"},
{"name":"pattern","label":"Pattern","type":"selection","values":[{"name":"TilePattern"},{"name":"ShiftTilePattern"}]},
{"name":"numDim","label":"Number of Dimensions","type":"range","value":2,"min":1,"max":5},
 {"type":"br"},
{"name":"dim_group","type":"group","quantity":"numDim","content":[
 {"name":"dim","label":"Extent Dimension","type":"logrange"},
 {"name":"dist","label":"Distribution","type":"selection",
  "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"TILE_blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"BLOCKCYCLIC_blockSize","label":"Block Size","type":"number","quantity":1}]}]}
 ]}
]
};

var ex02 = // how to not do it
{
"name": "ex02",
"type": "group",
"content": [
{"name":"numUnits","label":"Number of Units","type":"number"}, // Contained in Teamspec?
{"name":"pattern","label":"Pattern","type":"selection",
"values":[
{"name":"TilePattern<1>","params":[{"name":"size","label":"Size","type":"number","quantity":1},
                                   {"name":"dist","label":"Distribution","type":"selection","quantity":1,
                                    "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]}]},
                                   {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":1}
                                  ]},
{"name":"TilePattern<2>","params":[{"name":"size","label":"Size","type":"number","quantity":2},
                                   {"name":"dist","label":"Distribution","type":"selection","quantity":2,
                                    "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]}]},
                                   {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":2}
                                  ]},
{"name":"TilePattern<3>","params":[{"name":"size","label":"Size","type":"number","quantity":3},
                                   {"name":"dist","label":"Distribution","type":"selection","quantity":3,
                                    "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]}]},
                                   {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":3}
                                  ]},
{"name":"TilePattern<4>","params":[{"name":"size","label":"Size","type":"number","quantity":4},
                                   {"name":"dist","label":"Distribution","type":"selection","quantity":4,
                                    "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]}]},
                                   {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":4}
                                  ]},
{"name":"ShiftTilePattern<1>","params":[{"name":"size","label":"Size","type":"number","quantity":1},
                                        {"name":"dist","label":"Distribution","type":"selection","quantity":1,
                                         "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"blockSize","label":"Block Size","type":"number","quantity":1}]}]},
                                        {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":1}
                                  ]},
{"name":"ShiftTilePattern<2>","params":[{"name":"size","label":"Size","type":"number","quantity":2},
                                   {"name":"dist","label":"Distribution","type":"selection","quantity":2,
                                    "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"Block Size","type":"number","quantity":1}]}]},
                                   {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":2}
                                       ]},
{"name":"ShiftTilePattern<3>","params":[{"name":"size","label":"Size","type":"number","quantity":3},
                                        {"name":"dist","label":"Distribution","type":"selection","quantity":3,
                                         "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"Block Size","type":"number","quantity":1}]}]},
                                        {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":3}
                                  ]},
{"name":"ShiftTilePattern<4>","params":[{"name":"size","label":"Size","type":"number","quantity":4},
                                        {"name":"dist","label":"Distribution","type":"selection","quantity":4,
                                         "values":[{"name":"BLOCKED"},{"name":"CYCLIC"},{"name":"NONE"},{"name":"TILE","params":[{"name":"Block Size","type":"number","quantity":1}]},{"name":"BLOCKCYCLIC","params":[{"name":"Block Size","type":"number","quantity":1}]}]},
                                        {"name":"teamSpec","label":"TeamSpec","type":"number","quantity":4}
                                        ]}
]}
]
};
