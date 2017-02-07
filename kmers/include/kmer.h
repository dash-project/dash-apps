#ifndef KMER_
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

#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <string>

using ustring = std::basic_string<unsigned char>;

class KmerLine : public std::string {
  friend std::istream& operator>>(std::istream& is, KmerLine& line)
  {
    return std::getline(is, line);
  }
};

/* K-mer data structure */
typedef struct kmer {
  char sequence[KMER_PACKED_LENGTH];
  char l_ext;
  char r_ext;
} kmer_t;

typedef struct kmer_hash {
  std::size_t operator()(kmer_t const& kmer) const
  {
    return std::hash<std::string>{}(kmer.sequence);
  }
} kmer_hash_t;

bool operator==(const kmer_t& lhs, const kmer_t& rhs)
{
  return lhs.sequence == rhs.sequence && lhs.r_ext == rhs.r_ext &&
         lhs.l_ext == rhs.l_ext;
}

std::ostream& operator<<(std::ostream& stream, kmer_t const& kmer)
{
  stream << "seq: " << kmer.sequence << " r_ext: " << kmer.r_ext
         << " l_ext: " << kmer.l_ext;
  return stream;
}

#endif  // KMER_H
