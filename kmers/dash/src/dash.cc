#include <libdash.h>
#include <unistd.h>
#include <cassert>
#include <fstream>
#include <string>
#include <unordered_map>

// kmer specific files
#include <algos.h>
#include <kmer.h>

static dash::global_unit_t                                myid;
typedef dash::util::Timer<dash::util::TimeMeasure::Clock> Timer;

// verifies the file and returns the number of lines
size_t
verify_ufx_file(std::ifstream& ifs)
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
    std::cerr
        << "UFX text file has an unexpected line length for kmer length "
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

static std::string
read_line_nr(std::ifstream& ifs, size_t const nr)
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

static std::pair<int64_t, int64_t>
map_line_range(dash::global_unit_t const myid, size_t const nunits,
               size_t const nkmers)
{
  if (myid > nkmers) return std::make_pair(-1, -1);

  auto const rem        = nkmers % nunits;
  auto       nlocal     = nkmers / nunits;
  auto       line_start = myid * nlocal;

  if (rem) {
    if (myid < rem) {
      nlocal++;
      line_start = myid * nlocal;
    }
    else {
      auto const nbefore_rem = rem * (nlocal + 1);
      auto const nbefore     = (myid - rem) * nlocal;
      line_start             = nbefore_rem + nbefore;
    }
  }

  auto line_end = std::min(line_start + nlocal - 1, nkmers - 1);

  return std::make_pair(line_start, line_end);
}

static void
debug_lines(std::pair<int64_t, int64_t> const& lines, std::ifstream& file,
            dash::global_unit_t const myid)
{
  if (lines.first == -1 && lines.second == -1) return;

  std::string first_line = read_line_nr(file, lines.first);
  std::string last_line  = "";

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

int
main(int argc, char* argv[])
{
  typedef size_t              key_t;
  typedef kmer::packed_kmer_t mapped_t;
  typedef dash::UnorderedMap<key_t, mapped_t> map_t;
  typedef typename map_t::iterator   map_iterator;
  typedef typename map_t::value_type map_value;

  dash::init(&argc, &argv);
  Timer::Calibrate(0);

  auto const nunits = dash::size();
  myid              = dash::myid();

  dash::Shared<dash::default_size_t> nkmers{};
  dash::default_size_t nkmers_val = std::numeric_limits<std::size_t>::max();
  std::ifstream        ifs(argv[1]);

  if (0 == myid) {
    int wait = 1;
    while (wait)
      ;
    if (1 >= argc) {
      std::cout << "no input file given" << std::endl;
    }
    else {
      nkmers_val = verify_ufx_file(ifs);
    }
    nkmers.set(nkmers_val);
  }

  nkmers.barrier();
  nkmers_val = nkmers.get();

  // std::cout << "nlines: " << nkmers_val;

  if (nkmers_val == std::numeric_limits<std::size_t>::max()) {
    if (myid == 0) std::cout << "invalid input file!" << std::endl;
    dash::finalize();
    return EXIT_FAILURE;
  }

  auto line_range = map_line_range(myid, nunits, nkmers_val);

  if (line_range.first == -1) {
    dash::finalize();
    return EXIT_SUCCESS;
  }

  auto const nlocal =
      static_cast<map_t::size_type>(line_range.second - line_range.first + 1);

  map_t map{nkmers_val, nlocal / 4};

  map.barrier();

  // goto first local line
  ifs.seekg(line_range.first * LINE_SIZE, std::ios::beg);

  assert(ifs.good());

  /**********************************************************************
  1. Graph Construction
  ***********************************************************************/

  std::unordered_map<key_t, mapped_t> stdmap;
  Timer::timestamp_t ts_start = Timer::Now();

  transform_n(std::istream_iterator<kmer::KmerSequence>(ifs), nlocal,
              std::inserter(stdmap, stdmap.end()),
              [](kmer::KmerSequence const& input) {
                kmer::packed_kmer_t packed = pack_kmer(input);
                auto const          hasher = kmer::kmer_hash_t{};
                auto const          h      = hasher(packed);

                std::ostringstream os;
                os << "sequence: " << input.sequence() << ", hash: " << h
                   << ", lext: " << input.l_ext()
                   << ", r_ext: " << input.r_ext() << "\n";
                std::cout << os.str();
                return std::make_pair(h, packed);
              });

  Timer::timestamp_t duration = Timer::ElapsedSince(ts_start);

  std::cout << "Inserting " << nlocal
            << " Elements into a std::unordered_map took " << duration * 10e-3
            << " ms" << std::endl;

  ifs.clear();
  ifs.seekg(line_range.first * LINE_SIZE, std::ios::beg);

  assert(ifs.good());

  auto counter = 0;
  ts_start     = Timer::Now();
  transform_n(
      // iterator begin
      std::istream_iterator<kmer::KmerSequence>(ifs),
      // nelement
      nlocal,
      // insert iterator
      std::inserter(map, map.end()),
      // function
      [&counter, &nlocal](kmer::KmerSequence const& input) {
        // pack the kmer to decrease memory footprint
        kmer::packed_kmer_t packed = pack_kmer(input);
        // hash the packed kmer
        auto const h = kmer::kmer_hash_t{}(packed);

        //++counter;
        /*
        std::ostringstream os;
        os << "hash: " << h << ", lext: " << kmer.l_ext << ", r_ext: "
        << kmer.r_ext << "\n";
        std::cout << os.str();
        */
        //    std::cout << "---------------Inserted " << counter << "
        //    local elements" << "\n";
        return std::make_pair(h, packed);
      });

  duration = Timer::ElapsedSince(ts_start);

  std::cout << "Inserting " << nlocal
            << " Elements into a dash::UnorderedMap took " << duration * 10e-3
            << " ms" << std::endl;

  // 2. Filter all start nodes
  auto start_nodes =
      filter(map.lbegin(), map.lend(), [](map_t::value_type& pair) {
        //   std::cout << "Left ext: " << pair.second.l_ext << std::endl;
        return pair.second.l_ext == 'F';
      });

  //  std::cout <<  ": " << start_nodes.size() << "\n";

  /**********************************************************************
  2. Graph Traversal
  ***********************************************************************/

  if (myid == 0) {
    std::for_each(std::begin(start_nodes), std::end(start_nodes),
                  [](map_t::value_type const& pair) {
                    std::cout << "Start Node: " << pair.second.sequence
                              << "\n";
                  });
  }

  dash::finalize();
  return EXIT_SUCCESS;
}
