#ifndef PATTERN_VISUALIZER__input_param_h
#define PATTERN_VISUALIZER__input_param_h

#include <nlohmann/json.hpp>

#include <main.h>

using json = nlohmann::json;

typedef struct cli_params_t {
  cli_params_t() = default;
  bool output_blocks = true;
  bool output_memlayout = true;
  bool reduced = false;
} cli_params;

const cli_params default_output_params;

void print_basic_params();

void print_usage(char **argv);

//void print_params(const cli_params & params);

cli_params parse_args(int argc, char * argv[]);

json read_params();

const json getParamFromGroup(const json contents, const std::string name);

#endif // PATTERN_VISUALIZER__input_param_h
