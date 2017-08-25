#!/usr/bin/python

from datetime import datetime
import sys
import subprocess
import os
import math

import rpy2.robjects as robjects
robjects.r('library("scales")')  
import rpy2.robjects.lib.ggplot2 as ggplot2
ggplot2.theme_set(ggplot2.theme_bw ())
#print ggplot2.theme_get()
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector, StrVector, IntVector, DataFrame

def ggplot2_options ():
  return ggplot2.opts (**{'axis.title.x' : ggplot2.theme_blank(),
                          'axis.title.y' : ggplot2.theme_text(family = 'serif', face = 'bold', size = 15, angle=90, vjust=0.2),
                          'axis.text.x' : ggplot2.theme_text(family = 'serif', size = 15),
                          'axis.text.y' : ggplot2.theme_text(family = 'serif', size = 15),
                          'legend.title' : ggplot2.theme_text(family = 'serif', face = 'bold', size = 15),
                          'legend.text' : ggplot2.theme_text(family = 'serif', size = 15),
                          'aspect.ratio' : 0.6180339888,
    })

def ggplot2_colors ():
  return robjects.r('scale_fill_manual(values = c("#ffcb7e", "#1da06b", "#b94646", "#00368a"))')

def pdf_height (): return 3.7
def pdf_width (): return 7

pretty_varis = {"seq"      : "Sequential",
                "par"      : "Parallel",
                "expertseq": "Sequential (expert)",
                "expertpar": "Parallel (expert)"
                }
pretty_langs = {"chapel"   : "Chapel",
                "cilk"     : "Cilk",
                "erlang"   : "Erlang",
                "go"       : "Go",
                "scoop"    : "SCOOP",
                "tbb"      : "TBB"
                }

#languages = ["chapel", "cilk", "erlang", "go", "scoop", "tbb"]
languages = ["chapel", "cilk", "go", "tbb"]
problems = ["randmat", "thresh", "winnow", "outer", "product", "chain"]
variations = ["seq","par","expertseq","expertpar"]

start_actions = set(["start", "started", "restart", "resume"])
end_actions = set(["done", "pause"])

start_times = {}
total_times = {}


def main():
  #create_directories()
  read_table()
  load_data()
  #output_tables()
  #calculate()
  #test_significance()
  bargraph_variation()
  bargraph_variation_diff()
  bargraph_language()
  stat_test()
  simple_rank ()
  print_results ()

def print_results ():
  for lang in languages:
    for prob in problems:
      for var in variations:
        try:
          val = result[lang][prob][var]
        except KeyError:
          print (lang, prob, var)
          result[lang][prob][var] = 0.0
          
  for lang in languages:
    sys.stdout.write ("& " + pretty_langs [lang])
    for prob in ["randmat", "thresh", "winnow", "outer", "product", "chain"]:
      for var in ["seq", "expertseq", "par", "expertpar"]:
        if var.startswith('expert'):
          val = result[lang][prob][var] + result[lang][prob][var.replace('expert', '')]
        else:
          val = result[lang][prob][var]
        sys.stdout.write (" & " + str (int(round(val, 0))))
    #for var in ["seq", "expertseq", "par", "expertpar"]:
    #  sum = 0
    #  for prob in ["chain", "outer", "product", "randmat", "thresh", "winnow",]:
    #    if var.startswith('expert'):
    #      val = result[lang][prob][var] + result[lang][prob][var.replace('expert', '')]
    #    else:
    #      val = result[lang][prob][var]
    #    sum = sum + val
    #  sys.stdout.write (" & " + str (int(round(sum))))
    sys.stdout.write (" \\\\\n")

def simple_rank ():
  for lang in languages:
    for prob in problems:
      for var in variations:
        try:
          val = result[lang][prob][var]
        except KeyError:
          print (lang, prob, var)
          result[lang][prob][var] = 0.0

  for var in variations:
    print var
    for lang in languages:
      agg = 0
      for prob in problems:
        if var.startswith('expert'):
          val = result[lang][prob][var] + result[lang][prob][var.replace('expert', '')]
          valmin = min ([ result[l][prob][var] + result[l][prob][var.replace('expert', '')] for l in languages ])          
        else:
          val = result[lang][prob][var]
          valmin = min ([ result[l][prob][var] for l in languages ])
        agg = agg + float (val) / float (valmin)
      agg = agg / len (problems)
      print lang + '\t' + str (round (agg, 1))

def stat_test ():
  r = robjects.r

  res = {}

  for lang in languages:
    for prob in problems:
      for var in variations:
        try:
          val = result[lang][prob][var]
        except KeyError:
          print (lang, prob, var)
          result[lang][prob][var] = 0.0

  for var in variations:
    for lang1 in languages:
      if var.startswith('expert'):
        lang1_vals = FloatVector ([ result[lang1][prob][var] + result[lang1][prob][var.replace('expert', '')] for prob in problems ] )
      else:
        lang1_vals = FloatVector ([ result[lang1][prob][var] for prob in problems ] )
      for lang2 in languages:
        if var.startswith('expert'):
          lang2_vals = FloatVector ([ result[lang2][prob][var] + result[lang2][prob][var.replace('expert', '')] for prob in problems ] )
        else:
          lang2_vals = FloatVector ([ result[lang2][prob][var] for prob in problems ] )
        #print lang1
        #print lang1_vals
        #print lang2
        #print lang2_vals

        if lang1 not in res:
          res[lang1] = {}
        if lang2 not in res[lang1]:
          res[lang1][lang2] = {}

        pval = (r['wilcox.test'] (lang1_vals, lang2_vals, paired = True))[2][0]
        print pval
        if lang1 == lang2:
          res[lang1][lang2] = 0
        else:
          res[lang1][lang2] = pval

    # trivial data display
    print var
    print languages
    for lang1 in languages:
      sys.stdout.write (pretty_langs [lang1])
      for lang2 in languages:
        if lang1 == lang2:
          sys.stdout.write (" &        ")
        else:
          sys.stdout.write (" & " + str (round(res[lang1][lang2], 3)) + "  ")
      sys.stdout.write("\n")

def read_table():
  with open('tdist.txt', 'r') as f:
    linenum = 0
    for line in f:
      linenum += 1
      words = line.split()
      if linenum == 1:
        for i in range(len(words)):
          if i == 0:
            continue
          words[i] = words[i][0:-1]
          f = float(words[i])
          yindex.append(f)
      elif linenum == 2:
        continue;
      else:
        line = []
        for i in range(len(words)):
          if i == 0:
            xindex.append(float(words[i]))
          else:
            line.append(float(words[i]))
        table.append(line)
      #print words
  #print yindex
  #print table

def load_data():
  f = open("log_reverse.txt", "r")
  for line in f:
    words = line.split ()
    tz_offset = 5
    time_string = " ".join (words [:tz_offset])

    fmt = "%a %b %d %H:%M:%S %Y"
    parsed_date = datetime.strptime(time_string, fmt)
    #    print parsed_date
    words = words [tz_offset + 1:]
    #    print words
    if len(words) > 1:
      index = words[0]
      action = words[1]

      splits = index.split("-")

      if len(splits) > 1:
        language = splits[0]
        problem = splits[1]

        if language == "cpp" and len(splits) == 2 and problem != "refac":
          index = index + "-seq"
        elif problem == "chain":
          if len(splits) == 2:
            index = index + "-par"
      
      if len(index.split("-")) == 2 and index.split("-")[1] == "refac":
        pass
      elif action in start_actions:
        assert(index not in start_times);
        start_times[index] = parsed_date
      elif action in end_actions:
        end_time = parsed_date
        diff = end_time - start_times[index]
        del start_times[index]
        assert(diff.days == 0)
        diff = diff.seconds / 60.0
        if index in total_times:
          total_times[index] += diff
        else:
          total_times[index] = diff
      else:
        pass

  assert(len(start_times) == 0)

  for key, value in total_times.iteritems():
    words = key.split("-")
    if len(words) != 3:
      assert(False)
    language = words[0]
    problem = words[1]
    variation = words[2]

    if language not in result:
      result[language] = {}

    if problem not in result[language]:
      result[language][problem] = {}

    assert (variation not in result[language][problem])

    result[language][problem][variation] = value

  for problem in problems:
    if problem != "chain":
      result["tbb"][problem]["seq"] += result["cpp"][problem]["seq"]

  for language in languages:
    for problem in problems:
      if problem != "chain":
        result[language][problem]["par"] += result[language][problem]["seq"]

  # result["erlang"]["chain"]["seq"] = result["erlang"]["chain"]["par"]
  # result["scoop"]["chain"]["seq"] = result["scoop"]["chain"]["par"]

def bargraph_variation ():
  r = robjects.r
  for var in variations:
    # each variation gets plot
    values = []
    # normalized values
    nvalues = []

    langs = []
    probs = []

    for prob in problems:
      # aggregate by problems
      lvalues = []
      for lang in languages:
        # each problem displays a list of language times for that problem
 
        langs.append (pretty_langs [lang])
        probs.append (prob)
        value = 0
        try:
          value = result[lang][prob][var]
        except KeyError:
          print "Warning: no value for:"
          print (lang, prob, var)
          value = 0 # FIXME to account for missing seq-version of Erlang

        # for the expert times, add expert and non-expert times together
        if var.startswith('expert'):
          try:
            value = value + result[lang][prob][var.replace('expert','')]
          except KeyError:
            pass
        lvalues.append (value)
        
      values.extend (lvalues)
        
      lmin = min ([x for x in lvalues if x != 0])
      nvalues.extend ([(lambda x: x/lmin)(la) for la in lvalues])

    # plot histogram of actual times
    r.pdf ('bargraph-codingtime-var-' + var + '.pdf', height=pdf_height (), width=pdf_width ())

    df = robjects.DataFrame({'Language': StrVector (langs),
                             'Problem': StrVector (probs),
                             'Time' : FloatVector (values),
                             })

    dodge = ggplot2.position_dodge (width=0.9)
    gp = ggplot2.ggplot (df)

    pp = gp + \
        ggplot2.aes_string (x='Problem', y='Time', fill='Language') + \
        ggplot2.geom_bar (position='dodge', stat='identity') + \
        ggplot2_options () + \
        ggplot2_colors () + \
        robjects.r('scale_x_discrete(limits=c("randmat", "thresh", "winnow", "outer", "product", "chain"))') +\
        robjects.r('ylab("Coding time (in minutes)")')
 
    pp.plot ()

    # plot histogram of times normalized with respect to fastest time for a problem
    r.pdf ('bargraph-codingtime-var-norm-' + var + '.pdf', height=pdf_height (), width=pdf_width ())

    df = robjects.DataFrame({'Language': StrVector (langs),
                             'Problem': StrVector (probs),
                             'Time' : FloatVector (nvalues),
                             })

    dodge = ggplot2.position_dodge (width=0.9)
    gp = ggplot2.ggplot (df)

    pp = gp + \
        ggplot2.aes_string (x='Problem', y='Time', fill='Language') + \
        ggplot2.geom_bar (position='dodge', stat='identity') + \
        ggplot2_options () + \
        ggplot2_colors () + \
        robjects.r('scale_x_discrete(limits=c("randmat", "thresh", "winnow", "outer", "product", "chain"))') +\
        robjects.r('ylab("Coding time (normalized to fastest)")')

    pp.plot ()
    r['dev.off']()

def bargraph_variation_diff ():
  r = robjects.r

  for (standard, expert) in [('seq', 'expertseq'), ('par', 'expertpar')]:
    langs = []
    probs = []
    diffs  = []
    for lang in languages:
      for prob in problems:
        error = False
        try:
          time = result[lang][prob][standard]
        except KeyError:
          error = True
        try:
          time_expert = result[lang][prob][expert]
        except KeyError:
          error = True

        if not error:
          diff = (float(time_expert + time) / float(time) - 1)
        else:
          diff = 0

        langs.append (pretty_langs [lang])
        probs.append (prob)
        diffs.append (diff)

    r.pdf ('bargraph-codingtime-diff-' + standard + '.pdf', height=pdf_height (), width=pdf_width ())
    df = robjects.DataFrame({'Language': StrVector (langs),
                             'Problem': StrVector (probs),
                             'Difference' : FloatVector (diffs),
      })
    
    #print (df)
    gp = ggplot2.ggplot (df)
  
    pp = gp + \
        ggplot2.aes_string (x='Problem', y='Difference', fill='Language') + \
        ggplot2.geom_bar (position='dodge', stat='identity') + \
        ggplot2_options () + \
        ggplot2_colors () + \
        robjects.r('ylab("Coding time difference (in percent)")') +\
        robjects.r('scale_x_discrete(limits=c("randmat", "thresh", "winnow", "outer", "product", "chain"))') +\
        robjects.r('scale_y_continuous(labels = percent_format())')
    pp.plot ()
    r['dev.off']()

def bargraph_language ():
  r = robjects.r

  for language in languages:
    varis = []
    probs = []
    times  = []
    for prob in problems:
      for var in variations:
        try:
          time = result[language][prob][var]
        except KeyError:
          time = 0

        # for the expert times, add expert and non-expert times together
        if var.startswith('expert'):
          try:
            time = time + result[language][prob][var.replace('expert','')]
          except KeyError:
            pass
          
        varis.append (pretty_varis [var])
        probs.append (prob)
        times.append (time)
    r.pdf ('bargraph-codingtime-lang-' + language + '.pdf', height=pdf_height (), width=pdf_width ())
    df = robjects.DataFrame({'Variation': StrVector (varis),
                             'Problem': StrVector (probs),
                             'Time' : IntVector (times),
      })
    
    #print (df)
    gp = ggplot2.ggplot (df)
  
    pp = gp + \
        ggplot2.aes_string (x='Problem', y='Time', fill='Variation') + \
        ggplot2.geom_bar (position='dodge', stat='identity') + \
        ggplot2_options () + \
        ggplot2_colors () + \
        robjects.r('scale_x_discrete(limits=c("randmat", "thresh", "winnow", "outer", "product", "chain"))') +\
        robjects.r('ylab("Coding time (in minutes)")')
    pp.plot ()
    r['dev.off']()


    # Didn't change code below
    # ----------------------------------------------------------------------------------
    
def read_file_values(file_name):
  result = []
  f = open(file_name, 'r')
  for line in f:
    try:
      value = float(line)
      result.append(value)
    except ValueError:
      f.close()
      return []
  f.close()
  return result

def write_to_file(output, content):
  f = open(output, 'w')
  f.write(content)
  f.close()

def mean(x):
  n = len(x)
  total = 0
  for xi in x:
    total += xi
  return total / n

def stddev(x):
  n = len(x)
  total = 0
  meansq = mean([xi * xi for xi in x])
  meanx = mean(x)
  sqmean = meanx * meanx
  return math.sqrt((n * meansq - n * sqmean) / (n - 1))

yindex = []
xindex = []
table = []

def get_t(alpha, df):
  alpha *= 100
  y = -1
  for i in range(len(yindex)):
    if yindex[i] >= alpha:
      y = i
      break
  x = -1
  for i in range(len(xindex)):
    if xindex[i] > df:
      x = i
      break
  return table[x][y]

def ttest(xa, xb, alpha):
  '''meana = mean(xa)
  sa = stddev(xa)
  na = len(xa)
  delta = get_tdelta(xa, alpha)
  left = meana - delta
  right = meana + delta
  return left > 0'''
  conf = 1 - alpha
  with open('input.data', 'w') as f:
    f.write('x y\n')
    for i in range(len(xa)):
      f.write('%.10f\t%.10f\n' % (xa[i], xb[i]))
  cmd = '../r/test_paired.r %f > out.data' % conf
  # cmd = '../r/test_paired.r %f' % conf
  system(cmd)
  pvalue = read_file_values('out.data')[0]
  # print '%f ' % pvalue,
  # return pvalue <= alpha
  return pvalue

def get_tdelta(xa, alpha):
  meana = mean(xa)
  sa = stddev(xa)
  na = len(xa)
  t = get_t(1 - alpha / 2, na - 2)
  #print 'meana: %f\nsa2: %f\nna: %d\nt: %f' % (meana, sa * sa, na, t)
  return t * sa / math.sqrt(na)

result = {}
wc_result = {}
table_types = {"loc" : "-l", "now" : "-w"}

total_lines = 0

def system(cmd, timeout=False):
  ret = os.system(cmd)
  if ret != 0 and not timeout:
    print cmd
    assert(False)


def output_tables():
  def create_table(table_name, output_value, extra):
    old_stdout = sys.stdout

    for variation in variations:
      print variation
      sys.stdout = open(os.path.join (output_dir, "table-%s-%s.tex" % (
        table_name, variation)), "w")

      first = True
      print " & ",
      for problem in sorted(problems):
        if True:
          if not first:
            print " & ", 
          print problem,
          first = False;
      print " \\\\ \\hline"

      for language in sorted(languages):
        print language,
        for problem in sorted(problems):
          output_value(language, problem, variation, extra)
        print " \\\\"

    sys.stdout = old_stdout

  ########## time tables ###############

  def time_table_output(language, problem, variation, extra):
    print variation
    print (result[language][problem])
    assert(variation in result[language][problem])
    if variation in result[language][problem]:
      print " & ",
      print("%.2f" % result[language][problem][variation])

  create_table("time", time_table_output, None)

  ########## LoC-NoW-NoC tables ###############

  extensions = { "chapel" : "chpl", "cilk" : "cilk", "erlang" : "erl",
      "go" : "go", "scoop" : "e", "tbb" : "cc" }

  def wc_table_output(language, problem, variation, extra):
    extension = extensions[language]
    table_flag = extra["table_flag"]
    table_type = extra["table_type"]
    cmd = "find ../../%s/%s/%s/ | grep \"\\.%s$\" | xargs cat | grep . | wc %s > wc.out" % (
        language, problem, variation, extension, table_flag)
    if problem == "chain" and language in [
        "erlang"]:
      cmd = "find ../../%s/%s/ | grep \"\\.%s$\" | xargs cat | grep . | wc %s > wc.out" % (
          language, problem, extension, table_flag)
    os.system(cmd)

    if table_type not in wc_result:
      wc_result[table_type] = {}
    if language not in wc_result[table_type]:
      wc_result[table_type][language] = {}
    if problem not in wc_result[table_type][language]:
      wc_result[table_type][language][problem] = {}
    assert(variation not in wc_result[table_type][language][problem])
    value = open("wc.out", "r").read()

    wc_result[table_type][language][problem][variation] = value
    if table_type == "loc":
      global total_lines
      total_lines += int(value)
    print " & ", value,

  for table_name, table_flag in table_types.iteritems():
    create_table(table_name,  wc_table_output, {"table_flag": table_flag,
      "table_type" : table_name})

GRAPH_SIZE = 700
output_dir = os.path.abspath ("output/")
images_dir = os.path.join (output_dir, "images/")
chapters_dir = os.path.join (output_dir, "chapters/")

graph_dir = os.path.abspath ("graph/")


ALPHA = 0.1

ttest_res = {}

def output_pvalues(table, pretty, code):
  out = []
  out.append(
'''
\\begin{table}[htbp]
  \\caption{p-values for %s}
  \\label{tab:pv-%s}
  \\centering
  \\begin{tabular}{c|cccccc}
language''' % (pretty, code))
  for language in sorted(languages):
    out.append(' & %s' % language)
  out.append('\\\\\n\\hline\n')
  for la in sorted(languages):
    out.append('%s' % la)
    for lb in sorted(languages):
      if la == lb:
        out.append(' & --')
      else:
        out.append(' & %.3e' %  (table[la][lb]))
    out.append('\\\\\n')
  out.append('''
  \\end{tabular}
\\end{table}''')
  outstr = ''.join(out)
  output_file_name = "table-pvalue-%s.tex" % (code)
  output_file = os.path.join (chapters_dir, output_file_name)
  write_to_file(output_file, outstr)

def test_significance():
  types = ['tsing', 'tboth', 'size']
  for t in types:
    ttest_res[t] = {}

  for variation in variations:
    ttest_res['tsing'][variation] = {}
    for la in languages:
      ttest_res['tsing'][variation][la] = {}
      for lb in languages:
        if la == lb:
          continue
        left = []
        right = []
        for problem in problems:
          left.append(result[la][problem][variation])
          right.append(result[lb][problem][variation])
        pvalue = ttest(left, right, ALPHA)
        ttest_res['tsing'][variation][la][lb] = pvalue
        passed = pvalue <= ALPHA
        if passed:
          print '%s:%s:%s passed SINGLE' % (la, lb, variation)

  for la in languages:
    ttest_res['tboth'][la] = {}
    for lb in languages:
      if la == lb:
        continue
      left = []
      right = []
      for problem in problems:
        for variation in variations:
          left.append(result[la][problem][variation])
          right.append(result[lb][problem][variation])
      pvalue = ttest(left, right, ALPHA)
      ttest_res['tboth'][la][lb] = pvalue
      passed = pvalue <= ALPHA
      if passed:
        print '%s:%s passed BOTH' % (la, lb)

  for la in languages:
    ttest_res['size'][la] = {}
    for lb in languages:
      if la == lb:
        continue
      left = []
      right = []
      for problem in problems:
        for variation in variations:
          for table_type in table_types:
            left.append(int(wc_result[table_type][la][problem][variation]))
            right.append(
                int(wc_result[table_type][lb][problem][variation]))
      pvalue = ttest(left, right, ALPHA)
      ttest_res['size'][la][lb] = pvalue
      passed = pvalue <= ALPHA
      if passed:
        print '%s:%s passed SIZE' % (la, lb)

  output_pvalues(ttest_res['tsing']['seq'],
      'time to code sequential version', 'seq-ttc')
  output_pvalues(ttest_res['tsing']['par'],
      'time to code parallel version', 'par-ttc')
  output_pvalues(ttest_res['tboth'], 'time to code both versions',
      'both-ttc')
  output_pvalues(ttest_res['size'], 'size of source code', 'size')

def calculate():
  # left.append(result[la][problem][variation])
  variations = ['par']
  nmax = {}
  nmin = {}
  for problem in problems:
    nmax[problem] = {}
    nmin[problem] = {}
    for variation in variations:
      nmax[problem][variation] = 0
      nmin[problem][variation] = 999999999999
      for language in languages:
        nmax[problem][variation] = max(nmax[problem][variation],
            result[language][problem][variation])
        nmin[problem][variation] = min(nmin[problem][variation],
            result[language][problem][variation])
  calc = {}
  nsum = 0.0
  for language in languages:
    num = 0.0
    den = 0.0
    for problem in problems:
      for variation in variations:
        #num += result[language][problem][variation]
        #num -= nmin[problem][variation]
        #den += nmax[problem][variation]
        #den -= nmin[problem][variation]
        num += (result[language][problem][variation] / (
            float(nmin[problem][variation])))
        den += 1
    calc[language] = num / float(den)
    nsum += calc[language]

  vec = []
  for language in languages:
    #vec.append((calc[language] / nsum, language))
    vec.append((calc[language], language))

  print '\n\ntime to code (par)\n'

  vec.sort()
  for x in vec:
    (val, lang) = x
    print '%s : %f' % (lang, val)
  # int(wc_result[table_type][lb][problem][variation]))
  nmax = {}
  nmin = {}
  table_types = ['now']
  for t in table_types:
    nmax[t] = {}
    nmin[t] = {}
    for problem in problems:
      nmax[t][problem] = {}
      nmin[t][problem] = {}
      for variation in variations:
        nmax[t][problem][variation] = 0
        nmin[t][problem][variation] = 999999999999
        for language in languages:
          nmax[t][problem][variation] = max(nmax[t][problem][variation],
              int(wc_result[t][language][problem][variation]))
          nmin[t][problem][variation] = min(nmin[t][problem][variation],
              int(wc_result[t][language][problem][variation]))
  calc = {}
  nsum = 0.0
  for language in languages:
    num = 0.0
    den = 0.0
    for t in table_types:
      for problem in problems:
        for variation in variations:
          num += (int(wc_result[t][language][problem][variation]) / (
              float(nmin[t][problem][variation])))
          den += 1
    calc[language] = num / float(den)
    nsum += calc[language]

  vec = []
  for language in languages:
    #vec.append((calc[language] / nsum, language))
    vec.append((calc[language], language))

  print '\n\nsource code size (NoW - par)\n'

  vec.sort()
  for x in vec:
    (val, lang) = x
    print '%s : %f' % (lang, val)

def create_directories ():
  if not os.path.exists (output_dir):
    os.makedirs (output_dir)

  if not os.path.exists (chapters_dir):
    os.makedirs (chapters_dir)

  if not os.path.exists (images_dir):
    os.makedirs (images_dir)

if __name__ == "__main__":
  main()
