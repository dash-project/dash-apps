import unittest
import mox
import __builtin__
import main

class testMain(unittest.TestCase):
  def testGenerateErlangMain(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')
    main.write_to_file('directory/main.sh', 'dir=directory')

    m.StubOutWithMock(main, 'get_directory')
    main.get_directory('erlang', 'problem', 'variation').AndReturn(
        'directory')

    m.StubOutWithMock(main, 'get_problems_with_variations')
    main.get_problems_with_variations().AndReturn(
        [('problem', 'variation')])

    main.ERLANG_MAIN = 'dir=%s'

    m.ReplayAll()
    main.generate_erlang_main()
    m.UnsetStubs()
    m.VerifyAll()

  def testMakeAll(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'system')
    main.system('cd ../../language/problem/variation && make main')

    m.StubOutWithMock(main, 'get_all')
    main.get_all().AndReturn([('language', 'problem', 'variation')])

    m.ReplayAll()
    main.make_all()
    m.UnsetStubs()
    m.VerifyAll()

  def testCreateInputs(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')

    main.problem_classes = [
        main.RandmatProblem(), main.ThreshProblem(), main.WinnowProblem(),
        main.OuterProblem(), main.ProductProblem(), main.ChainProblem()]
    main.inputs = [main.ProblemInput(10, 15, 20, 30, 40)]

    main.write_to_file('randmat_10_15_20.in', '10 15 20\n')
    main.write_to_file('thresh_10_15_30.in', '10 15 30\n')
    main.write_to_file('winnow_10_15_30_40.in', '10 15 40\n')
    main.write_to_file('outer_40.in', '40\n')
    main.write_to_file('product_40.in', '40\n')
    main.write_to_file('chain_10_20_30_40.in', '10\n20\n30\n40\n')

    m.ReplayAll()
    main.create_inputs()
    m.UnsetStubs()
    m.VerifyAll()

  def testRunAll(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'system')
    m.StubOutWithMock(main, 'read_from_file')
    m.StubOutWithMock(main, 'get_all')

    m.StubOutWithMock(main, 'get_directory')
    main.get_directory('language', 'problem', 'variation').AndReturn(
        'directory')

    main.get_all().AndReturn([('language', 'problem', 'variation')])

    main.inputs = [main.ProblemInput(10, 15, 20, 30, 40)]
    main.TIMEOUT = 99
    main.threads = [999]
    main.system(('' # 'timeout 99 ' # TODO timeout on ensisun
                 'taskset -c 0-998 '
                 '/usr/bin/time -a -f %e -o '
                 'time-language-problem-variation-0-999.out '
                 'directory/main < '
                 'problem_10_15_20_30_40.in > /dev/null 1>&0 2>&0'),
                 timeout=True)
    main.read_from_file('time-language-problem-variation-0-999.out')

    m.ReplayAll()
    main.run_all()
    m.UnsetStubs()
    m.VerifyAll()

  def testGetResults(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'get_all')
    main.get_all().AndReturn([('language', 'problem', 'variation')])

    m.StubOutWithMock(main, 'read_file_values')
    main.read_file_values(
        'time-language-problem-variation-0-999.out').AndReturn([1, 2])

    main.threads = [999]

    m.ReplayAll()
    main.get_results()
    self.assertEqual(
        main.results[999]['problem']['variation']['language'][0],
        (1 + 2) / 2.)
    m.UnsetStubs()
    m.VerifyAll()

  def testOutputExecTimeGraphs(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')
    m.StubOutWithMock(main, 'system')

    main.problems = ['problem']
    main.variations = ['seq']
    main.languages = ['language']

    main.threads = [999]
    main.results = {
        999: {'problem': {
            'seq': {'language': {0: 100}}}}}

    main.write_to_file('output/images/graph-exec-time-seq-0.perf', '=cluster;language\ncolors=black,yellow,red,med_blue,light_green,cyan\n=table\nyformat=%g\n=norotate\nxscale=1\nmax=150.000000\nylabel=Sequential execution time in seconds for input 0\nproblem 100.0000000000\n')
    main.system('output/bargraph.pl -fig output/images/graph-exec-time-seq-0.perf | fig2dev -L ppm -m 4 > output/images/graph-exec-time-seq-0.ppm')
    main.system('mogrify -reverse -flatten output/images/graph-exec-time-seq-0.ppm')
    main.system('mogrify -resize 700x700 -format png output/images/graph-exec-time-seq-0.ppm')
    main.write_to_file('output/chapters/graph-exec-time-0.tex',
        '\\begin{figure}[htbp]\n  %\\centering\n  \\includegraphics[width=125mm]{images/graph-exec-time-seq-0.png}\n  \\caption{Sequential Execution Time}\n  \\label{fig:exec:time:seq:0}\n\\end{figure}\n')

    m.ReplayAll()
    main.create_graph("exec-time", main.results[999], '', use_subfigure=False)
    m.UnsetStubs()
    m.VerifyAll()

  def testOutputSpeedupGraphs(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')
    m.StubOutWithMock(main, 'system')

    main.problems = ['problem']
    main.variations = ['seq', 'par']
    main.languages = ['language']

    main.threads = [2, 4]
    main.results = {
        2: {'problem': {
            'par': {'language': {0: 50}}}},
        4: {'problem': {
            'seq': {'language': {0: 100}},
            'par': {'language': {0: 25}}}}}

    main.write_to_file(
        'output/images/graph-speedup-language-problem-0.dat', '2\t50.0000000000\t2.0000000000\t2\t1.0000000000\t1\n4\t25.0000000000\t4.0000000000\t4\t1.0000000000\t1\n')

    main.system('cp output/images/graph-speedup-language-problem-0.dat plot.dat')
    main.system('gnuplot output/plot.script')
    main.system('mv plot.png output/images/graph-speedup-language-problem-0.png')
    main.system('rm plot.dat')
    main.write_to_file('output/chapters/graph-speedup-language-problem-0.tex',
        '\\begin{figure}[htbp]\n  %\\centering\n  \\includegraphics[width=100mm]{images/graph-speedup-language-problem-0.png}\n  \\caption{Speedup and efficiency for language language in problem problem}\n  \\label{fig:exec:spd:language:problem:0}\n\\end{figure}\n')
    main.write_to_file('output/chapters/graph-speedup.tex', '\\input{chapters/graph-speedup-language-problem-0.tex}\n')

    m.ReplayAll()
    main.create_speedup_graph("speedup", main.results, use_subfigure=False)
    m.UnsetStubs()
    m.VerifyAll()

  def testOutputProblemSpeedupGraphs(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')
    m.StubOutWithMock(main, 'system')

    main.problems = ['problem']
    main.variations = ['seq', 'par']
    main.languages = ['language']

    main.threads = [2, 4]
    main.results = {
        2: {'problem': {
            'par': {'language': {0: 50}}}},
        4: {'problem': {
            'seq': {'language': {0: 100}},
            'par': {'language': {0: 25}}}}}

    main.write_to_file('other.script',   '\nset xrange [0:8]\nset xtics 1\nset yrange [0:8]\nset ytics 1\nset xlabel "threads"\nset ylabel "speedup"\nset terminal png\nset output "plot.png"\nset key left\nplot \'output/images/graph-speedup-language-problem-0.dat\' using 1:4 title \'ideal speedup\' w lp, \'output/images/graph-speedup-language-problem-0.dat\' using 1:3 title \'language speedup\' w lp')
    main.system('gnuplot other.script')
    main.system('mv plot.png output/images/graph-problem-speedup-problem-0.png')
    main.system('rm other.script')
    main.write_to_file('output/chapters/graph-problem-speedup-problem-0.tex',  '\\begin{figure}[htbp]\n  %\\centering\n  \\includegraphics[width=100mm]{images/graph-problem-speedup-problem-0.png}\n  \\caption{Speedup for problem problem in all languages}\n  \\label{fig:exec:spd:problem:0}\n\\end{figure}\n')
    main.write_to_file('output/chapters/graph-problem-speedup.tex',  '\\input{chapters/graph-problem-speedup-problem-0.tex}\n')

    m.ReplayAll()
    main.create_problem_speedup_graph("problem-speedup", "speedup", use_subfigure=False)
    m.UnsetStubs()
    m.VerifyAll()

  def testOutputLanguageSpeedupGraphs(self):
    m = mox.Mox()

    m.StubOutWithMock(main, 'write_to_file')
    m.StubOutWithMock(main, 'system')

    main.problems = ['problem']
    main.variations = ['seq', 'par']
    main.languages = ['language']

    main.threads = [2, 4]
    main.results = {
        2: {'problem': {
            'par': {'language': {0: 50}}}},
        4: {'problem': {
            'seq': {'language': {0: 100}},
            'par': {'language': {0: 25}}}}}

    main.write_to_file('other.script',   '\nset xrange [0:8]\nset xtics 1\nset yrange [0:8]\nset ytics 1\nset xlabel "threads"\nset ylabel "speedup"\nset terminal png\nset output "plot.png"\nset key left\nplot \'output/images/graph-speedup-language-problem-0.dat\' using 1:4 title \'ideal speedup\' w lp, \'output/images/graph-speedup-language-problem-0.dat\' using 1:3 title \'problem speedup\' w lp')
    main.system('gnuplot other.script')
    main.system('mv plot.png output/images/graph-language-speedup-language-0.png')
    main.system('rm other.script')
    main.write_to_file('output/chapters/graph-language-speedup-language-0.tex',  '\\begin{figure}[htbp]\n  %\\centering\n  \\includegraphics[width=100mm]{images/graph-language-speedup-language-0.png}\n  \\caption{Speedup for language language in all problems}\n  \\label{fig:exec:spd:language:0}\n\\end{figure}\n')
    main.write_to_file('output/chapters/graph-language-speedup.tex',  '\\input{chapters/graph-language-speedup-language-0.tex}\n')

    m.ReplayAll()
    main.create_language_speedup_graph("language-speedup", "speedup", use_subfigure=False)
    m.UnsetStubs()
    m.VerifyAll()

if __name__ == '__main__':
  unittest.main()
