#ifndef EXTRAE_INSTRUMENTATION__H
#define EXTRAE_INSTRUMENTATION__H


#ifdef USE_EXTRAE
#include "extrae_user_events.h"
#include "extrae_types.h"
extern "C" void Extrae_init (void) __attribute__((weak));
extern "C" void Extrae_event (extrae_type_t type, extrae_value_t value) __attribute__((weak));
extern "C" void Extrae_fini (void) __attribute__((weak));
extern "C" void Extrae_define_event_type (extrae_type_t *type, char *type_description, unsigned *nvalues, extrae_value_t *values, char **values_description) __attribute__((weak));
#endif

#define EVENT_NONE  0
#define EVENT_POTRF 1
#define EVENT_TRSM  2
#define EVENT_GEMM  3
#define EVENT_SYRK  4
#define EVENT_PREFETCH 5

#ifdef USE_EXTRAE
static extrae_type_t et;
static extrae_value_t ev[6] = {0, 10, 11, 12, 13, 20};
static char *extrae_names[] = {"none", "potrf()", "trsm()", "gemm()", "syrk()", "PREFETCH"};

#define EXTRAE_ENTER(_e) Extrae_event(et, ev[_e])
#define EXTRAE_EXIT(_e)  Extrae_event(et, ev[EVENT_NONE])
#else
#define EXTRAE_ENTER(_e)
#define EXTRAE_EXIT(_e)
#endif

#endif // EXTRAE_INSTRUMENTATION__H
