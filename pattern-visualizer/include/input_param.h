#ifndef PATTERN_VISUALIZER__input_param_h
#define PATTERN_VISUALIZER__input_param_h

#include <main.h>

typedef struct cli_params_t {
  cli_params_t()
  { }
  std::string type    = "summa";
  extent_t    size_x  = 110;
  extent_t    size_y  = 110;
  extent_t    units_x = 10;
  extent_t    units_y = 10;
  int         tile_x  = -1;
  int         tile_y  = -1;
  bool        blocked_display = false;
  bool        balance_extents = false;
  bool        cout    = false;
} cli_params;

static const
cli_params default_params;

void print_usage(char **argv);

void print_params(const cli_params & params);

cli_params parse_args(int argc, char * argv[]);

#endif // PATTERN_VISUALIZER__input_param_h
