#ifndef LOGGING_H
#define LOGGING_H

#ifdef ENABLE_LOGGING

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>

#include <dash/Init.h>
#include <dash/Team.h>

#define LOG(msg) (Log(__FILE__, __LINE__, LogData<None>() << msg))

// Workaround GCC 4.7.2 not recognizing noinline attribute
#ifndef NOINLINE_ATTRIBUTE
#ifdef __ICC
#define NOINLINE_ATTRIBUTE __attribute__((noinline))
#else
#define NOINLINE_ATTRIBUTE
#endif  // __ICC
#endif  // NOINLINE_ATTRIBUTE

#define PRINT_MYID(os)                          \
  do {                                          \
    os << "[UNIT " << std::setw(4) << dash::myid() << "] "; \
  } while (0)

struct None {
};

template <typename List>
struct LogData {
  List list;
};

template <typename List>
void Log(const char* file, int line, LogData<List>&& data) NOINLINE_ATTRIBUTE
{
  std::ostringstream os;
  PRINT_MYID(os);
  os << "[" << std::setw(6) << file << ":" << line << "] ";
  output(os, std::move(data.list));
  os << std::endl;
  std::cout << os.str();
}

template <typename Begin, typename Value>
constexpr LogData<std::pair<Begin&&, Value&&>> operator<<(
    LogData<Begin>&& begin, Value&& value) noexcept
{
  return {{std::forward<Begin>(begin.list), std::forward<Value>(value)}};
}

template <typename Begin, size_t n>
constexpr LogData<std::pair<Begin&&, const char*>> operator<<(
    LogData<Begin>&& begin, const char (&value)[n]) noexcept
{
  return {{std::forward<Begin>(begin.list), value}};
}

typedef std::ostream& (*PfnManipulator)(std::ostream&);

template <typename Begin>
constexpr LogData<std::pair<Begin&&, PfnManipulator>> operator<<(
    LogData<Begin>&& begin, PfnManipulator value) noexcept
{
  return {{std::forward<Begin>(begin.list), value}};
}

template <typename Begin, typename Last>
void output(std::ostream& os, std::pair<Begin, Last>&& data)
{
  output(os, std::move(data.first));
  os << data.second;
}

inline void output(std::ostream& os, None)
{
}

#define MAX_ELEMS_RANGE_LOGGING 100

#define LOG_TRACE_RANGE(desc, begin, end)                                   \
  do {                                                                      \
    using value_t =                                                         \
        typename std::iterator_traits<decltype(begin)>::value_type;         \
    using difference_t =                                                    \
        typename std::iterator_traits<decltype(begin)>::difference_type;    \
    auto const nelems = std::distance(begin, end);                          \
    auto const max_elems =                                                  \
        std::min<difference_t>(nelems, MAX_ELEMS_RANGE_LOGGING);            \
    std::ostringstream os;                                                  \
    os << desc << " : ";                                                    \
    std::copy(                                                              \
        begin, begin + max_elems, std::ostream_iterator<value_t>(os, " ")); \
    if (nelems > MAX_ELEMS_RANGE_LOGGING) os << "...";                      \
    LOG(os.str());                                                          \
  } while (0)

#else

#define LOG(msg) \
  do {           \
  } while (0)

#define LOG_TRACE_RANGE(desc, begin, end) \
  do {                                    \
  } while (0)

#endif

#endif
