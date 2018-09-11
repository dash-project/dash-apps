#ifndef EXTRAE__H
#define EXTRAE__H

#define USE_EXTRAE


enum ExtraeEventsEnum {
  EVENT_NONE = 0,
  TIME_INCREMENT,
  CALCFORCEFORNODES,
  SYNCFORCE,
  CALCACCELERATIONFORNODES,
  CALCVELPOSFORNODES,
  LAGRANGEELEMS,
  MATERIALPROPERTIES,
  MONOQGRADIENTS,
  MONOQELEMS,
  COURANTCONSTRAINTS,
  HYDROCONSTRAINTS,
  _NUM_EVENTS
};

enum ExtraeTaskEventsEnum {
  EVENT_TASK_NONE = 0,
  HOURGLASSCONTROLFORELEM,
  FBHOURGLASSFORCEFORELEMS,
  FBHOURGLASSFORCEFORELEMSREDUCTION,
  INITSTRESSTERMSFORELEMS,
  INTEGRATESTRESSFORELEMS,
  INTEGRATESTRESSFORELEMSREDUCTION,
  CHECKDETERM,
  RELEASE,
  ALLOCATESTRAINS,
  ALLOCATEGRADIENTS,
  _NUM_TASK_EVENTS
};

//#include "extrae_user_events.h"
//#include "extrae_types.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_EXTRAE

typedef unsigned extrae_type_t;
typedef unsigned long long extrae_value_t;
void Extrae_init (void) __attribute__((weak));
void Extrae_event (extrae_type_t type, extrae_value_t value) __attribute__((weak));
void Extrae_define_event_type (extrae_type_t *type, const char *type_description, unsigned *nvalues, extrae_value_t *values, const char **values_description) __attribute__((weak));
void Extrae_fini (void) __attribute__((weak));
static extrae_type_t et = 1000;
static extrae_value_t extrae_ids[_NUM_EVENTS] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110};
static const char *extrae_names[_NUM_EVENTS] = {
  "NONE", "TimeIncrement", "CalcForceForNodes", "SyncForce",
  "CalcAccelerationForNodes", "CalcVelPosForNodes", "LagrangeElems",
  "MaterialProperties", "CalcMonotonicQGradientsForElems",
  "CalcMonotonicQForElems", "CourantConstraintForElems",
  "HydroConstraintForElems"};

static extrae_type_t ett = 1001;
static extrae_value_t extrae_task_ids[_NUM_TASK_EVENTS] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
static const char *extrae_names_tasks[_NUM_TASK_EVENTS] = {
  "NONE", "CalcHourglassControlForElems", "CalcFBHourglassForceForElems",
  "CalcFBHourglassForceForElems_REDUCTION", "InitStressTermsForElems",
  "IntegrateStressForElems", "IntegrateStressForElems_REDUCTION", "CheckDeterm",
  "RELEASE", "AllocateStrains", "AllocateGradients"};

#ifdef __cplusplus
}
#endif

#define EXTRAE_ENTER(_e) if (Extrae_event) Extrae_event(et, extrae_ids[_e])
#define EXTRAE_EXIT(_e)  if (Extrae_event) Extrae_event(et, extrae_ids[EVENT_NONE])
#define EXTRAE_INIT() do { \
  if (Extrae_define_event_type) { \
    unsigned nvalues = _NUM_EVENTS; \
    Extrae_define_event_type(&et, "Lulesh Functions", &nvalues, extrae_ids, extrae_names); \
  } \
} while (0)
#define EXTRAE_TASK_ENTER(_e) if (Extrae_event) Extrae_event(ett, extrae_task_ids[_e])
#define EXTRAE_TASK_EXIT(_e)  if (Extrae_event) Extrae_event(ett, extrae_task_ids[EVENT_NONE])
#define EXTRAE_TASK_INIT() do { \
  if (Extrae_define_event_type) { \
    unsigned nvalues = _NUM_TASK_EVENTS; \
    Extrae_define_event_type(&ett, "Lulesh Task Functions", &nvalues, extrae_task_ids, extrae_names_tasks); \
  } \
} while (0)

#else  // USE_EXTRAE
#define EXTRAE_ENTER(_e)
#define EXTRAE_EXIT(_e)
#define EXTRAE_INIT()
#endif // USE_EXTRAE

struct extrae_event {
  extrae_event(ExtraeEventsEnum ev) : ev{ev} {
    EXTRAE_ENTER(ev);
  }
  ~extrae_event(){
    EXTRAE_EXIT(ev);
  }
  ExtraeEventsEnum ev;
};

struct extrae_task_event {
  extrae_task_event(ExtraeTaskEventsEnum ev) : ev{ev} {
    EXTRAE_TASK_ENTER(ev);
  }
  ~extrae_task_event(){
    EXTRAE_TASK_EXIT(ev);
  }
  ExtraeTaskEventsEnum ev;
};


#endif // EXTRAE__H


