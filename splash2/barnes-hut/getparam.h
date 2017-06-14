#ifndef _GETPARAM_H_
#define _GETPARAM_H_

#include "stdinc.h"

#ifdef __cplusplus
extern "C" {
#endif

void initparam(const char ** defv);
string getparam(string name);
long getiparam(string name);
long getlparam(string name);
bool getbparam(string name);
double getdparam(string name);
long scanbind(const char**, string name);
bool matchname(string, string name);
string extrvalue(string arg);

#ifdef __cplusplus
}
#endif

#endif
