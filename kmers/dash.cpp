#include <contig_generation.h>
#include <libdash.h>
#include <fstream>
#include <string>
#include <cassert>

// verifies the file and returns the number of lines
size_t verify_ufx_file(std::ifstream& ifs)
{
  size_t err_val = std::numeric_limits<std::size_t>::max();
  if (ifs.bad()) return err_val;

  // save state
  std::istream::iostate state_backup = ifs.rdstate();
  // clear state
  ifs.clear();
  // Pos Backup
  std::istream::streampos pos_backup = ifs.tellg();

  ifs.seekg(std::ios::beg);

  std::string line;
  std::getline(ifs, line);

  if (line.size() != (LINE_SIZE - 1)) {
    std::cerr << "UFX text file has an unexpected line length for kmer length "
              << line.size();
    return err_val;
  }

  if (line.at(KMER_LENGTH) != ' ' && line.at(KMER_LENGTH) != '\t') {
    std::cerr << "Unexpected format for firstLine: " << line;
    return err_val;
  }

  ifs.seekg(std::ios::beg);

  size_t nlines = std::count(std::istreambuf_iterator<char>(ifs),
                             std::istreambuf_iterator<char>(), '\n');

  // if the file is not end with '\n' , then line count  should plus 1
  ifs.unget();
  if (ifs.get() != '\n') {
    ++nlines;
  }
  // recover state
  ifs.clear();  // previous reading may set eofbit
  ifs.seekg(pos_backup);
  ifs.setstate(state_backup);

  return nlines;
}
int main(int argc, char* argv[])
{
  typedef int key_t;
  typedef double mapped_t;
  typedef dash::UnorderedMap<key_t, mapped_t> map_t;
  typedef typename map_t::iterator map_iterator;
  typedef typename map_t::value_type map_value;

  dash::init(&argc, &argv);

  auto const nunits = dash::size();
  auto const myid = dash::myid();

  dash::Shared<dash::default_size_t> nkmers{};
  dash::default_size_t nkmers_val = std::numeric_limits<std::size_t>::max();
  std::ifstream file(argv[1]);

  if (0 == myid) {
    if (1 >= argc) {
      std::cout << "no input file given" << std::endl;
    }
    else {
      nkmers_val = verify_ufx_file(file);
    }
    nkmers.set(nkmers_val);
  }

  nkmers.barrier();
  nkmers_val = nkmers.get();

  if (nkmers_val == std::numeric_limits<std::size_t>::max()) {
    if (myid == 0) std::cout << "invalid input file!" << std::endl;
    dash::finalize();
    return EXIT_FAILURE;
  }

  if (myid == 0) std::cout << " we have " << nkmers_val << " lines in our file" << std::endl;

  map_t map{nkmers_val, 1};

  map.barrier();

  auto const lcap = map.lcapacity();
  //size_t line_start = myid * lcap;
  size_t line_start = 10;
  size_t line_end = line_start + lcap;

  //goto specific line
  file.seekg((line_start - 1) * LINE_SIZE, std::ios::beg);

  assert(file.good());

  std::string kmer;

  std::getline(file, kmer);

  std::cout << "Line " << line_start << " is: " << kmer << std::endl;


  dash::finalize();
  return EXIT_SUCCESS;
}