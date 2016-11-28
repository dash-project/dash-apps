#ifndef LOGGER_H
#define LOGGER_H

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>

#define LOG_UNIT(logger) \
     BOOST_LOG_TRIVIAL(logger)
#endif // LOGGER_H
