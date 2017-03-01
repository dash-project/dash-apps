var ex01 = // todo complete
{
"name": "ex01",
"type": "group",
"content": [
{"name":"pattern","label":"Pattern","type":"selection","values":[{"name":"summa","label":"SummaPattern"},{"name":"block","label":"BlockPattern"},{"name":"shift","label":"ShiftTilePattern"},{"name":"seq","label":"SeqTilePattern"}]},
{"name":"numDim","label":"Number of Dimensions","type":"range","value":2,"min":2,"max":2},
 {"type":"br"},
{"name":"blocked_display","label":"Blocked display","type":"checkbox"},
{"name":"balance_extents","label":"Balance extents","type":"checkbox"},
 {"type":"br"},
{"name":"dim_group","type":"group","quantity":"numDim","content":[
 {"name":"size","label":"Extent Dimension","type":"logrange"},
 {"name":"units","label":"Units Dimension","type":"number"},
 {"name":"tile","label":"Tile Dimension","type":"number"}
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
