#line 185 "/home/kowalewski/workspaces/splash2/codes/apps/barnes/../../null_macros/c.m4.null.POSIX_BARRIER"

#line 1 "grav.H"
#ifndef _GRAV_H_
#define _GRAV_H_

void hackgrav(bodyptr p, long ProcessId);
void gravsub(nodeptr const p, long ProcessId);
void hackwalk(long ProcessId);
void walksub(nodeptr n, real dsq, long ProcessId);
bool subdivp(nodeptr const p, real dsq, long ProcessId);

#endif
