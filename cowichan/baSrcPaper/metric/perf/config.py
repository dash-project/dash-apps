from problems import *

languages = ["chapel", "cilk", "go", "tbb"]
#languages = ["scoop"]
problems = ["randmat", "thresh", "winnow", "outer", "product", "chain"]
#problems = ["randmat", "thresh", "winnow", "outer", "product"]
#variations = ["seq","expertseq","expertpar"]
variations = ["seq","expertseq","expertpar", "par"]
threads = [1, 2, 4, 8, 16, 32]
output_dir = "output"


inputs = [
    #ProblemInput(20000, 20000, 666, 1, 1),
    #ProblemInput(20000, 20000, 666, 1, 10000),
    ProblemInput(32, 8000, 666, 1, 10000),
  ]

class Config:
  def __init__ (self):
    self.languages = languages
    self.problems = problems
    self.variations = variations
    self.inputs = inputs
    self.threads = threads
    self.output_dir = output_dir

cfg = Config ()

def get_problems_with_variations():
  for problem in sorted(cfg.problems):
    for variation in sorted(cfg.variations):
      yield (problem, variation)

def get_all():
  for (problem, variation) in get_problems_with_variations():
    for language in sorted(cfg.languages):
      yield (language, problem, variation)
