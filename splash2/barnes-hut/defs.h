/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

#ifndef _DEFS_H_
#define _DEFS_H_

//#include <assert.h>
#include <dash/Atomic.h>
#include <dash/GlobRef.h>
#include "stdinc.h"
#include "vectmath.h"

#define MAX_PROC 128
#define MAX_BODIES_PER_LEAF 10
/* maximum number of locks: Originally 2048 but MPI does crash with such a
 * high number */
#define MAXLOCK 128    /* maximum number of locks on DASH */
#define PAGE_SIZE 4096 /* in bytes */

#define NSUB (1 << NDIM) /* subcells per cell */

/*
 * BODY and CELL data structures are used to represent the tree:
 *
 *         +-----------------------------------------------------------+
 * root--> | CELL: mass, pos, cost, quad, /, o, /, /, /, /, o, /, done |
 *         +---------------------------------|--------------|----------+
 *                                           |              |
 *    +--------------------------------------+              |
 *    |                                                     |
 *    |    +--------------------------------------+         |
 *    +--> | BODY: mass, pos, cost, vel, acc, phi |         |
 *         +--------------------------------------+         |
 *                                                          |
 *    +-----------------------------------------------------+
 *    |
 *    |    +-----------------------------------------------------------+
 *    +--> | CELL: mass, pos, cost, quad, o, /, /, o, /, /, o, /, done |
 *         +------------------------------|--------|--------|----------+
 *                                       etc      etc      etc
 */

/*
 * NODE: data common to BODY and CELL structures.
 */

typedef dash::GlobPtr<struct _node> nodeptr;
typedef struct _node *nodelptr;

typedef struct _node {
  // 8
  long type; /* code for node type: body or cell */
  // 8
  real mass; /* total mass of node */
  // 24
  vector pos; /* position of node */
  // 8
  long cost; /* number of interactions computed */
  // 8
  long level;
  // 40
  nodeptr parent; /* ptr to parent of this node in tree */
  // 8
  long child_num; /* Index that this node should be put at in parent cell */
} node;

/*
#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Cost(x) (((nodeptr) (x))->cost)
#define Level(x) (((nodeptr) (x))->level)
#define Parent(x) (((nodeptr) (x))->parent)
#define ChildNum(x) (((nodeptr) (x))->child_num)
*/

#define Type(x) ((x).member<long>(&node::type))
#define Mass(x) ((x).member<real>(&node::mass))
#define Pos(x) ((x).member<vector>(&node::pos))
#define Cost(x) ((x).member<long>(&node::cost))
#define Level(x) ((x).member<long>(&node::level))
#define Parent(x) ((x).member<nodeptr>(&node::parent))
#define ChildNum(x) ((x).member<long>(&node::child_num))

/*
 * BODY: data structure used to represent particles.
 */

typedef dash::GlobPtr<struct _body> bodyptr;
typedef struct _body *bodylptr;
typedef dash::GlobPtr<struct _leaf> leafptr;
typedef dash::GlobPtr<struct _cell> cellptr;

static_assert(sizeof(nodeptr) == sizeof(cellptr),
              "GlobPtr of two types do not equal in sizeof");

#define BODY static_cast<long>(01) /* type code for bodies */

struct _body_old {
  //-------- Common to both body and cell
  long type;
  real mass;  /* mass of body */
  vector pos; /* position of body */
  long cost;  /* number of interactions computed */
  long level;
  leafptr parent;
  long child_num; /* Index that this node should be put */

  //-------- Unique attribute in body
  vector vel; /* velocity of body */
  vector acc; /* acceleration of body */
  real phi;   /* potential at body */
};

// 56
typedef struct _body : public node {
  //-------- Unique attribute in body
  // 24
  vector vel; /* velocity of body */
  // 24
  vector acc; /* acceleration of body */
  // 8
  real phi; /* potential at body */
} body;

std::ostream &operator<<(std::ostream &os, const body &b);

static_assert(sizeof(struct _body_old) == sizeof(struct _body),
              "struct _body does not satisfy standard_layout");

/*
#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)
*/
#define Vel(x) (x.member<vector>(&body::vel))
#define Acc(x) (x.member<vector>(&body::acc))
#define Phi(x) (x.member<real>(&body::phi))

/*
 * CELL: structure used to represent internal nodes of tree.
 */

#define CELL static_cast<long>(02) /* type code for cells */

struct _cell_old {
  //-------- Common to both body and cell
  long type;
  real mass;      /* total mass of cell */
  real pos[NDIM]; /* cm. position of cell */
  long cost;      /* number of interactions computed */
  long level;
  cellptr parent;
  long child_num; /* Index [0..8] that this node should be put */
                  //-------- Unique attribute in cell
  long processor; /* Used by partition code */
  // Not really used!!! So no problem in DASH
  struct _cell *next, *prev; /* Used in the partition array */
  long seqnum;
#ifdef QUADPOLE
  matrix quad; /* quad. moment of cell */
#endif
  /*volatile*/ long done; /* flag to tell when the c.of.m is ready */
  nodeptr subp[NSUB];     /* descendents of cell */
};

typedef struct _cell : public node {
  //-------- Unique attribute in cell
  long processor; /* Used by partition code */
  // Not really used!!! So no problem in DASH
  struct _cell *next, *prev; /* Used in the partition array */
  long seqnum; // seqnum for locks
#ifdef QUADPOLE
  matrix quad; /* quad. moment of cell */
#endif
  /*volatile*/ long done; /* flag to tell when the c.of.m is ready */
  nodeptr subp[NSUB];     /* descendents of cell */
} cell;

std::ostream &operator<<(std::ostream &os, const cell &c);

static_assert(sizeof(struct _cell) == sizeof(struct _cell_old),
              "struct _cell does not satisfy standard_layout");

//#define Subp(x) (((cellptr) (x))->subp)
#define Subp(x) ((x).member<nodeptr>(&cell::subp))
#define SeqnumC(x) ((x).member<long>(&cell::seqnum))

/*
 * LEAF: structure used to represent leaf nodes of tree.
 */

#define LEAF static_cast<long>(03) /* type code for leaves */

struct _leaf_old {
  cellptr parent;
  long type;
  real mass;  /* total mass of leaf */
  vector pos; /* cm. position of leaf */
  long cost;  /* number of interactions computed */
  long level;
  long child_num; /* Index [0..8] that this node should be put */

  long processor; /* Used by partition code */
  // Not really used!!! So no problem in DASH
  struct _leaf_old *next, *prev; /* Used in the partition array */
  long seqnum;
#ifdef QUADPOLE
  matrix quad; /* quad. moment of leaf */
#endif
  volatile long done; /* flag to tell when the c.of.m is ready */
  long num_bodies;
  bodyptr bodyp[MAX_BODIES_PER_LEAF]; /* bodies of leaf */
};

typedef struct _leaf : public node {
  long processor; /* Used by partition code */
  // Not really used!!! So no problem in DASH
  struct _leaf *next, *prev; /* Used in the partition array */
  long seqnum;
#ifdef QUADPOLE
  matrix quad; /* quad. moment of leaf */
#endif
  /*volatile */ long done; /* flag to tell when the c.of.m is ready */
  long num_bodies;
  bodyptr bodyp[MAX_BODIES_PER_LEAF]; /* bodies of leaf */
} leaf;

static_assert(sizeof(struct _leaf_old) == sizeof(struct _leaf),
              "struct _leaf does not satisfy standard_layout");

#define SeqnumL(x) ((x).member<long>(&leaf::seqnum))

namespace std {
template <>
class is_standard_layout<leaf> : public integral_constant<bool, true> {
};
template <>
class is_standard_layout<cell> : public integral_constant<bool, true> {
};
template <>
class is_standard_layout<body> : public integral_constant<bool, true> {
};
}

//#define Bodyp(x)  (((leafptr) (x))->bodyp)
#define Bodyp(x) (x.member<nodeptr>(&cell::subp))

#ifdef QUADPOLE
// TODO:
#define Quad(x) (((cellptr)(x))->quad)
#endif
#define Done(x) ((x).member<long>(&cell::done))
#define AtomicDone(x) dash::GlobRef<dash::Atomic<long>>(Done(x).dart_gptr())

/*
typedef dash::GlobRef<dash::Atomic<long>> atomic_flag_t;
typedef dash::GlobRef<long> flag_t;

template<typename T>
inline dash::GlobRef<dash::Atomic<T>> Atomic(dash::GlobRef<T> flag) {
  return dash::GlobRef<dash::Atomic<T>>(flag.dart_gptr());
}
*/

/* SOME TYPE CONVERSION */
#define NODE_AS_CELL(x) (static_cast<cellptr>(x))
#define NODE_AS_LEAF(x) (static_cast<leafptr>(x))

/*
 * Integerized coordinates: used to mantain body-tree.
 */
#define MAXLEVEL ((8L * (long)sizeof(long)) - 2L)
#define IMAX (1L << MAXLEVEL) /* highest bit of int coord */
#endif
