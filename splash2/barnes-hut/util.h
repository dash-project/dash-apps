#ifndef _UTIL_H_
#define _UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

double xrand(double xl, double xh);
void pranset(long seed);
double prand(void);
double cputime(void);

#ifdef __cplusplus
}
#endif

#endif
