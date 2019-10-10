
#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <time.h>

#define TIMESTAMP(time_) 						\
  do {									\
      struct timespec ts;						\
      clock_gettime(CLOCK_MONOTONIC, &ts);				\
      time_=((double)ts.tv_sec)+(1.0e-9)*((double)ts.tv_nsec);		\
  } while(0)


#endif /* UTIL_H_INCLUDED */

