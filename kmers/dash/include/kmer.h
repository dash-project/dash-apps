#ifndef KMER_H
#define KMER_H

#ifndef MAXIMUM_CONTIG_SIZE
#define MAXIMUM_CONTIG_SIZE 100000
#endif

#ifndef KMER_LENGTH
#define KMER_LENGTH 19
#endif

#ifndef LOAD_FACTOR
#define LOAD_FACTOR 1
#endif

#ifndef LINE_SIZE
#define LINE_SIZE (KMER_LENGTH + 4)
#endif

#ifndef KMER_PACKED_LENGTH
#define KMER_PACKED_LENGTH 5
#endif

#include <algorithm>
#include <cassert>
#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

using ustring = std::basic_string<unsigned char>;

namespace kmer {

class KmerSequence {
  friend std::istream& operator>>(std::istream& is, KmerSequence& line);

  public:
    KmerSequence(std::string & sequence, char l_ext, char r_ext)
    : sequence_(sequence),
      l_ext_(l_ext),
      r_ext_(r_ext)
    {}

    KmerSequence() = default;
    KmerSequence(KmerSequence const & other) = default;
    KmerSequence(KmerSequence && other) = default;
    KmerSequence & operator=(KmerSequence const & other) = default;
    KmerSequence & operator=(KmerSequence && other) = default;
    ~KmerSequence() = default;

  public:
    std::string const sequence() const { return sequence_; }
    char l_ext() const { return l_ext_; }
    char r_ext() const { return r_ext_; }

  private:
    std::string sequence_;
    char l_ext_;
    char r_ext_;
};

std::istream& operator>>(std::istream& is, KmerSequence& line)
{
  std::string content_;
  auto& ret = std::getline(is, content_);

  if (!content_.empty()) {
    assert(content_.size() == LINE_SIZE - 1);
    line.sequence_ = content_.substr(0, KMER_LENGTH);
    line.l_ext_ = content_[KMER_LENGTH + 1];
    line.r_ext_ = content_[KMER_LENGTH + 2];
  }

  return ret;
}

/* K-mer data structure */
typedef struct packed_kmer {
  char sequence[KMER_PACKED_LENGTH];
  char l_ext;
  char r_ext;

  packed_kmer() = default;
  packed_kmer(char const * seq, char l_ext, char r_ext)
    : l_ext(l_ext),
      r_ext(r_ext)
  {
    std::memcpy(sequence, seq, sizeof(KMER_PACKED_LENGTH));
  }
} packed_kmer_t;

typedef struct kmer_hash {
  std::size_t operator()(packed_kmer_t const& kmer) const
  {
    size_t hashval = 5381;
    for (int i = 0; i < KMER_PACKED_LENGTH; i++) {
      hashval = kmer.sequence[i] + (hashval << 5) + hashval;
    }
    return hashval;
  }
} kmer_hash_t;

static unsigned char convertFourMerToPackedCode(const unsigned char* fourMer)
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

static std::string const pack_sequence(std::string const& seq_to_pack)
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

packed_kmer_t pack_kmer(kmer::KmerSequence const & kmer) {
  auto packed_seq = pack_sequence(kmer.sequence());
  packed_kmer_t packed{packed_seq.c_str(), kmer.l_ext(), kmer.r_ext()};
  return packed;
}

bool operator==(const packed_kmer_t& lhs, const packed_kmer_t& rhs)
{
  return std::memcmp(lhs.sequence, rhs.sequence, KMER_PACKED_LENGTH) == 0 &&
         lhs.r_ext == rhs.r_ext && lhs.l_ext == rhs.l_ext;
}

std::ostream& operator<<(std::ostream& stream, packed_kmer_t const& kmer)
{
  stream << "seq: " << kmer.sequence << " r_ext: " << kmer.r_ext
         << " l_ext: " << kmer.l_ext;
  return stream;
}

}


#endif  // KMER_H
