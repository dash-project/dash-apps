#ifndef KMER_ALGOS_H
#define KMER_ALGOS_H

#include <algorithm>

template <typename InputIterator, typename Size, typename OutputIterator,
          typename UnaryOperation>
OutputIterator transform_n(InputIterator first, Size n, OutputIterator result,
                           UnaryOperation op)
{
  return std::generate_n(result, n, [&first, &op]() { return op(*first++); });
}

template <class InputIterator, class Functor>
auto filter(InputIterator begin, InputIterator end, Functor f)
    -> std::vector<typename std::iterator_traits<InputIterator>::value_type>
{
  using ValueType = typename std::iterator_traits<InputIterator>::value_type;

  std::vector<ValueType> result;
  result.reserve(std::distance(begin, end));

  std::copy_if(begin, end, std::back_inserter(result), f);

  return result;
}
#endif
