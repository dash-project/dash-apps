#!/usr/bin/env python
import os
import math
import time
from subprocess import Popen, call, PIPE

from problems import *
from utils import *
from config import *

ERLANG_MAIN = ("#!/bin/sh\n"
               "cd ~/lucia/metric/perf/%s\n"
               "erl -noshell +S $1 -s main main is_bench -s init stop\n")

def generate_erlang_main():
  for (problem, variation) in get_problems_with_variations():
    directory = get_directory("erlang", problem, variation)
    output = directory + "/main.sh"
    print output
    write_to_file(output, ERLANG_MAIN % directory)
    cmd = "chmod +x %s" % output
    system(cmd)

def make_all():
  for (language, problem, variation) in get_all():
    directory = get_directory(language, problem, variation)
    cmd = "cd %s && make main" % directory
    print cmd
    system(cmd)

def create_inputs():
  for i in range(len(problem_classes)):
    problem = problem_classes[i]
    for j in range(len(inputs)):
      cur = inputs[j]
      file_name = problem.input_file_name(cur)
      write_to_file(file_name, problem.get_input(cur))

def run_all(redirect_output=True):
  # TODO: check processor usage
  for nthreads in threads:
    for (language, problem, variation) in get_all():
      if is_sequential (variation) and nthreads != threads[-1]: continue
      for i in range(len(inputs)):
        #TODO: get time output file name
        time_output = get_time_output(language, problem, variation, i, nthreads)
        # TODO: refactor variations
        #print time_output
        cmd = ""
        env = ""
        if language == "go" and is_parallel (variation):
          env += "GOMAXPROCS=%d " % nthreads

        cmd += "taskset "
        if is_sequential (variation) or nthreads == 1:
          cmd += "-c 0 "
        else:
          cmd += "-c 0-%d " % (nthreads - 1)

        # using python timing below
        # directory = get_directory(language, problem, variation)
        # cmd += "/usr/bin/time -a -f %%e -o %s %s/" % (time_output,
        #     directory)
          
        directory = get_directory(language, problem, variation)
        # cmd += directory + '/'

        if language == "erlang":
          cmd += 'erl -noshell -pa ' + directory + '/ -s main main is_bench -s init stop'
        elif language == "scoop":
          cmd += directory + '/' + "main -bench -i"
        else:
          cmd += directory + '/' + "main"

        if language == "chapel":
          cmd += " --numLocales=1"
          if is_parallel (variation):
            cmd += " --numThreadsPerLocale=%d" % (nthreads)
          else:
            cmd += " --numThreadsPerLocale=1"
        elif language == "cilk":
          if is_parallel (variation):
            cmd += " --nproc %d" % nthreads
          else:
            pass # must NOT pass --nproc here (because of --is_bench)
        elif language == 'tbb':
          if is_parallel (variation):
            cmd += ' --threads %d' % nthreads
        elif language == 'erlang':
          if is_parallel (variation):
            cmd += ' +S %d' % nthreads
          else:
            cmd += ' +S 1'

        if (language == "chapel" or language == "cilk" or language == "tbb"
            or language == 'go'):
          cmd += " --is_bench"

        #if language != "scoop":
        #  cmd += " <";

        # cmd += " %s" % (problem_map[problem].input_file_name(inputs[i]));

        #if redirect_output:
        #  cmd += " > /dev/null 1>&0 2>&0"

        cmd = env + '/usr/bin/time -f "%M" ' + cmd
        print cmd

        t1 = time.time ()
        with open (problem_map[problem].input_file_name(inputs[i]), 'r') as input_file:
          proc = Popen (cmd, stdin=input_file, stdout=PIPE, stderr=PIPE, shell=True)
          (out, mem_usage) = proc.communicate ()
        t2 = time.time ()
        tdiff = t2 - t1

        # divide GNU time memory measurement by 4 if it's GNU time V1.7

        with open (get_mem_output (language, problem, variation), 'w') as output_file:
          output_file.write (str (float (mem_usage.strip()) / 4))
        with open (time_output, 'a') as output_file:
          output_file.write (str (tdiff) + '\n')

        print tdiff

TOTAL_EXECUTIONS = 2

def main():
  print "Building all programs"
  make_all()
  print "Creating inputs"
  create_inputs()
  print "Running tests"
  for _ in range(TOTAL_EXECUTIONS):
    run_all(redirect_output=False)  # TODO: remove outputs

if __name__ == '__main__':
  main()
