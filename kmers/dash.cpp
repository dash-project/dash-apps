#include <unistd.h>
#include <libdash.h>
#include <fstream>
#include <string>
#include <cassert>
#include <unordered_map>

#include <kmer.h>

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

static std::string read_line_nr(std::ifstream& ifs, size_t const nr)
{
  // save state
  std::istream::iostate state_backup = ifs.rdstate();
  // clear state
  ifs.clear();
  // Pos Backup
  std::istream::streampos pos_backup = ifs.tellg();
  // goto specific line
  ifs.seekg(nr * LINE_SIZE, std::ios::beg);

  assert(ifs.good());

  std::string ret;
  std::getline(ifs, ret);

  // recover state
  ifs.clear();  // previous reading may set eofbit
  ifs.seekg(pos_backup);
  ifs.setstate(state_backup);

  return ret;
}

static std::pair<int64_t, int64_t> map_line_range(
    dash::global_unit_t const myid, size_t const nunits, size_t const nkmers)
{
  if (myid > nkmers) return std::make_pair(-1, -1);

  auto const rem = nkmers % nunits;
  auto nlocal = nkmers / nunits;
  auto line_start = myid * nlocal;

  if (rem) {
    if (myid < rem) {
      nlocal++;
      line_start = myid * nlocal;
    } else {
      auto const nbefore_rem = rem * (nlocal + 1);
      auto const nbefore = (myid - rem) * nlocal;
      line_start = nbefore_rem + nbefore;
    }
  }

  auto line_end = std::min(line_start + nlocal - 1, nkmers - 1);

  return std::make_pair(line_start, line_end);
}

static void debug_lines(std::pair<int64_t, int64_t> const& lines,
                        std::ifstream& file, dash::global_unit_t const myid)
{
  if (lines.first == -1 && lines.second == -1) return;

  std::string first_line = read_line_nr(file, lines.first);
  std::string last_line = "";

  if (lines.first < lines.second) {
    last_line = read_line_nr(file, lines.second);
  }

  std::stringstream ss;

  if (first_line != "")
    ss << "[Unit " << std::setw(2) << myid << "]: First Line " << std::setw(2)
       << (lines.first + 1) << " is: " << first_line << "\n";
  if (last_line != "") {
    ss << "[Unit " << std::setw(2) << myid << "]: Last Line  " << std::setw(2)
       << (lines.second + 1) << " is: " << last_line << "\n";
  }

  std::cout << ss.str();
}

template <typename InputIterator, typename OutputIterator,
          typename UnaryOperation>
OutputIterator transform_n(InputIterator first, size_t n, OutputIterator result,
                           UnaryOperation op)
{
  return std::generate_n(result, n, [&first, &op]() { return op(*first++); });
}

unsigned char convertFourMerToPackedCode(const unsigned char* fourMer)
{
  int retval = 0;
  int code, i;
  int pow = 64;

  for (i = 0; i < 4; i++) {
    char base = fourMer[i];
    switch (base) {
      case 'A':
        code = 0;
        break;
      case 'C':
        code = 1;
        break;
      case 'G':
        code = 2;
        break;
      case 'T':
        code = 3;
        break;
    }
    retval += code * pow;
    pow /= 4;
  }
  return ((unsigned char)retval);
}

std::string const pack_sequence(std::string const& seq_to_pack)
{
  int ind, j = 0;  // coordinate along unpacked string ( matches with m_len )
  int i = 0;       // coordinate along packed string

  size_t m_len = seq_to_pack.size();
  assert(m_len == KMER_LENGTH);

  ustring const us(seq_to_pack.begin(), seq_to_pack.end());

  ustring m_data;

  // do the leading seq in blocks of 4
  for (; j <= m_len - 4; i++, j += 4) {
    auto const sub = us.substr(j, 4);
    auto const converted = convertFourMerToPackedCode(sub.c_str());
    m_data.push_back(converted);
  }

  // last block is special case if m_len % 4 != 0: append "A"s as filler
  int const remainder = m_len % 4;
  std::string const tmp = "AAAA";
  ustring blockSeq(tmp.begin(), tmp.end());

  for (ind = 0; ind < remainder; ind++) {
    blockSeq.at(ind) = us.at(j + ind);
  }

  m_data.push_back(convertFourMerToPackedCode(blockSeq.c_str()));
  assert(m_data.size() == KMER_PACKED_LENGTH);

  return std::string(m_data.begin(), m_data.end());
}

int main(int argc, char* argv[])
{
  typedef int key_t;
  typedef kmer_t mapped_t;
  typedef dash::UnorderedMap<key_t, mapped_t> map_t;
  typedef typename map_t::iterator map_iterator;
  typedef typename map_t::value_type map_value;

  dash::init(&argc, &argv);

  auto const nunits = dash::size();
  auto const myid = dash::myid();

  dash::Shared<dash::default_size_t> nkmers{};
  dash::default_size_t nkmers_val = std::numeric_limits<std::size_t>::max();
  std::ifstream ifs(argv[1]);

  if (0 == myid) {
    std::cout << "my id is: " << ::getpid() << std::endl;
    if (1 >= argc) {
      std::cout << "no input file given" << std::endl;
    } else {
      nkmers_val = verify_ufx_file(ifs);
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

  map_t map{nkmers_val, 1};

  map.barrier();

  auto line_range = map_line_range(myid, nunits, nkmers_val);

  if (line_range.first == -1) {
    dash::finalize();
    return EXIT_SUCCESS;
  }

  auto const nlocal = line_range.second - line_range.first + 1;
  std::vector<std::string> lines(nlocal);

  // goto first local line
  ifs.seekg(line_range.first, std::ios::beg);

  assert(ifs.good());

  std::unordered_map<std::string, kmer_t> map_kmer;

  /**********************************************************************
  1. Graph Construction
  ***********************************************************************/

  // 1. Map all kmers from file to unordered map
  transform_n(std::istream_iterator<KmerLine>(ifs), nlocal,
              std::inserter(map_kmer, std::end(map_kmer)),
              [](KmerLine const & input) {

    std::cout << "length of line: " << input.size();
    assert(input.size() == LINE_SIZE - 1);

    std::string const seq = input.substr(0, KMER_LENGTH);
    char const l_ext = input[KMER_LENGTH+1];
    char const r_ext = input[KMER_LENGTH+2];
    std::string const packed = pack_sequence(seq);

    kmer_t kmer{};
    packed.copy(kmer.sequence, KMER_PACKED_LENGTH);
    kmer.l_ext = l_ext;
    kmer.r_ext = r_ext;

    return std::make_pair(packed, kmer);
  });

  // 2. Filter all start nodes
  std::vector<std::pair<std::string, kmer_t>> start_nodes;
  start_nodes.reserve(std::distance(std::begin(map_kmer), std::end(map_kmer)));

  std::copy_if(std::begin(map_kmer), std::end(map_kmer),
               std::back_inserter(start_nodes),
               [](std::pair<std::string, kmer_t> const & item) {
    return item.second.l_ext == 'F';
  });


  /**********************************************************************
  2. Graph Traversal
  ***********************************************************************/

  if (myid == 0) {
    std::for_each(map_kmer.begin(), map_kmer.end(),
                  [](std::pair<std::string, kmer_t> const& iter) {
      std::cout << "Original: " << iter.second.sequence
                << ", packed: " << iter.first << "\n";
    });
    std::for_each(start_nodes.begin(), start_nodes.end(),
                  [](std::pair<std::string, kmer_t> const& iter) {
      std::cout << "Start Node: " << iter.second.sequence << "\n";
    });
  }

  dash::finalize();
  return EXIT_SUCCESS;
}