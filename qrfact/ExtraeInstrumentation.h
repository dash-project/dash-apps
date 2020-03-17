#ifndef EXTRAE_INSTRUMENTATION__H
#define EXTRAE_INSTRUMENTATION__H


#ifdef USE_EXTRAE
#include "extrae_user_events.h"
#include "extrae_types.h"
extern "C" void Extrae_init (void);
extern "C" void Extrae_event (extrae_type_t type, extrae_value_t value);
extern "C" void Extrae_fini (void);
extern "C" void Extrae_define_event_type (extrae_type_t *type, char *type_description, unsigned *nvalues, extrae_value_t *values, char **values_description);
#endif

#define EVENT_NONE  0
#define EVENT_GEQRT 1
#define EVENT_ORMQR 2
#define EVENT_TSQRT 3
#define EVENT_TSMQR 4
#define EVENT_WRKK  5
#define EVENT_WRKN  6

#ifdef USE_EXTRAE
static extrae_type_t et;
static extrae_value_t ev[7] = {0, 10, 11, 12, 13, 14, 15};
static char *extrae_names[7] = {"none", "GEQRT", "ORMQR", "TSQRT", "TSMQR", "WRKK", "WRKN"};

#define EXTRAE_ENTER(_e) Extrae_event(et, ev[_e])
#define EXTRAE_EXIT(_e)  Extrae_event(et, ev[EVENT_NONE])
#else
#define EXTRAE_ENTER(_e)
#define EXTRAE_EXIT(_e)
#endif

#endif // EXTRAE_INSTRUMENTATION__H
