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

#include <pthread.h>

#include <sys/time.h>

#include <unistd.h>

#include <stdlib.h>

#include <malloc.h>

#include "Logging.h"
#include "code.h"
#include "load.h"

extern pthread_t PThreadTable[];

std::ostream &operator<<(std::ostream &os, const body &b)
{
  os << "{"
     << "type: " << b.type << ", "
     << "mass: " << std::setprecision(7) << std::fixed << b.mass << ", "
     << "pos: [" << std::setprecision(5) << std::fixed << b.pos[0] << ", "
     << b.pos[1] << ", " << b.pos[2] << "], "
     << "cost: " << b.cost << ", "
     << "level: " << b.level
     << ", "
     //<< "parent: " << b.parent << ", "
     << "child_num: " << b.child_num << ", "
     << "vel: [" << b.vel[0] << ", " << b.vel[1] << ", " << b.vel[2] << "], "
     << "acc: [" << b.acc[0] << ", " << b.acc[1] << ", " << b.acc[2] << "], "
     << "phi: " << b.phi << "}";

  return os;
};

std::ostream &operator<<(std::ostream &os, const cell &c)
{
  os << "{"
     << "type: " << c.type << ", "
     << "mass: " << std::setprecision(7) << c.mass << ", "
     << "pos: [" << std::scientific << c.pos[0] << ", " << c.pos[1] << ", "
     << c.pos[2] << "], "
     << "cost: " << c.cost << ", "
     << "level: " << c.level << ", "
     << "child_num: " << c.child_num << ", "
     << "parent: " << c.parent << ", "
     << "processor: " << c.processor << ", "
     << "seq_num: " << c.seqnum << ", "
     << "subp: [";

  for (auto idx = 0; idx < NSUB; ++idx) {
    os << c.subp[idx] << ",";
  }

  os << "]}";

  return os;
}

std::ostream &operator<<(std::ostream &os, const leaf &c)
{
  os << "{"
     << "type: " << c.type << ", "
     << "mass: " << std::setprecision(7) << c.mass << ", "
     << "pos: [" << std::scientific << c.pos[0] << ", " << c.pos[1] << ", "
     << c.pos[2] << "], "
     << "cost: " << c.cost << ", "
     << "level: " << c.level << ", "
     << "child_num: " << c.child_num << ", "
     << "parent: " << c.parent << ", "
     << "processor: " << c.processor << ", "
     << "seq_num: " << c.seqnum << ", "
     << "num_bodies: " << c.num_bodies << "}";
#if 0
     << "subp: [";

  for (auto idx = 0; idx < NSUB; ++idx) {
    os << c.subp[idx] << ",";
  }

  os << "]}";
#endif

  return os;
}

/*
 * MAKETREE: initialize tree structure for hack force calculation.
 */

void maketree(long ProcessId)
{
  Local.myncell = 0;
  Local.mynleaf = 0;
  if (ProcessId == 0) {
    // Process 0 initially gets the root
    Local.mycelltab[Local.myncell++] = G_root.get();
  }
  // The processes remotely access the root
  Local.Current_Root = static_cast<cellptr>(G_root.get());

  for (auto pp = std::begin(Local.mybodytab);
       pp < std::next(std::begin(Local.mybodytab), Local.mynbody);
       ++pp) {
    auto       bodyp = *pp;
    auto const pmass = Mass(*bodyp).get();
    // insert all bodies into the tree
    if (pmass != 0.0) {
      Local.Current_Root = loadtree(
          bodyp, static_cast<cellptr>(Local.Current_Root), ProcessId);
    }
    else {
      /*
      {pthread_mutex_lock(&(Global->io_lock));};v
      fprintf(stderr, "Process %ld found body %ld to have zero mass\n",
              ProcessId, (long)p);
      {pthread_mutex_unlock(&(Global->io_lock));};
      */
      ASSERT(false);
    }
  }
  dash::barrier();

  LOG("mynleaf: " << Local.mynleaf << ", myncell: " << Local.myncell);

  hackcofm(ProcessId);

  dash::barrier();
}

cellptr InitCell(cellptr parent, long ProcessId)
{
  auto c = makecell(ProcessId);
  ASSERT(c.local());
  auto &c_val     = *(c.local());
  c_val.processor = ProcessId;
  c_val.next      = NULL;
  c_val.prev      = NULL;
  if (parent == cellptr{nullptr}) {
    c_val.level = IMAX >> 1;
  }
  else {
    c_val.level = Level((*parent)) >> 1;
  }

  c_val.parent    = parent;
  c_val.child_num = 0;
  return c;
}

leafptr InitLeaf(cellptr parent, long ProcessId)
{
  // allocate a local leaf
  auto leaf_gptr = makeleaf(ProcessId);
  ASSERT(leaf_gptr.local());
  auto &lv     = *(leaf_gptr.local());
  lv.processor = ProcessId;
  lv.next      = NULL;
  lv.prev      = NULL;
  if (parent == cellptr{nullptr}) {
    lv.level = IMAX >> 1;
  }
  else {
    lv.level = Level((*parent)) >> 1;
  }
  lv.parent = parent;

  lv.child_num = 0;

  return (leaf_gptr);
}

static void _printtree(nodeptr n, std::ostringstream &ss)
{
  long    k;
  cellptr c;
  cell    c_val;
  leaf    l_val;
  leafptr l;
  bodyptr p;
  // long nseq;
  long const type = Type(*n);
  switch (type) {
    case CELL:
      c     = NODE_AS_CELL(n);
      c_val = static_cast<cell>(*c);
      ss << "Cell : Cost = " << c_val.cost << ", ";
      PRTV(ss, "Pos", c_val.pos);
      ss << "\n";
      for (k = 0; k < NSUB; k++) {
        ss << "Child #" << k << ": ";
        if (c_val.subp[k] == nodeptr{}) {
          ss << "NONE";
        }
        else {
          if (Type(*(c_val.subp[k])) == CELL) {
            ss << "C: Cost = " << static_cast<long>(Cost(*(c_val.subp[k])));
          }
          else {
            auto const leaf_val = (*NODE_AS_LEAF(c_val.subp[k])).get();
            ss << "L: # Bodies = " << leaf_val.num_bodies
               << ", Cost = " << leaf_val.cost;
          }
          auto tmp = (*(c_val.subp[k])).get();
          ss << ", ";
          PRTV(ss, "Pos", tmp.pos);
        }
        ss << "\n";
      }
      for (k = 0; k < NSUB; k++) {
        if (c_val.subp[k] != nodeptr{}) {
          _printtree(c_val.subp[k], ss);
        }
      }
      break;
    case LEAF:
      l     = NODE_AS_LEAF(n);
      l_val = static_cast<leaf>(*l);
      ss << "Leaf : # Bodies = " << l_val.num_bodies
         << ", Cost = " << l_val.cost << ", ";
      PRTV(ss, "Pos", l_val.pos);
      ss << "\n";
      for (k = 0; k < l_val.num_bodies; k++) {
        p = l_val.bodyp[k];

        // BUG: We cannot copy around dash::GlobPtr and use all functionality.
        // This is because a dash::GlobPtr has some native pointer members
        // with invalid addresses on a different unit. For this reason we
        // construct a new global pointer with the proper members based on a
        // DART pointer
        //auto const p_tmp = bodyptr{bodytab.globmem(), p.dart_gptr()};

        auto const p_val = (*p).get();
        ss << "Body #" << p - static_cast<bodyptr>(bodytab.begin())
           << ": Num = " << k << ", Level = " << p_val.level << ", ";
        PRTV(ss, "Pos", p_val.pos);
        ss << "\n";
      }
      break;
    default:
      std::cerr << "Bad type\n";
      exit(-1);
      break;
  }
}

void printtree(nodeptr n)
{
  std::ostringstream ss;
  _printtree(n, ss);
  ss << std::endl;
  std::cout << ss.str();
}

/*
 * LOADTREE: descend tree and insert particle.
 */

nodeptr loadtree(bodyptr p, cellptr root, long ProcessId)
{
  long    l, xp[NDIM], _xor[NDIM];
  bool    flag;
  long    i, j, root_level;
  bool    valid_root;
  long    kidIndex;
  cellptr mycell;
  leafptr le;

  auto body_val = static_cast<body>(*p);
  auto root_val = static_cast<cell>(*root);
  // get integerized coords of real position in space
  intcoord(xp, body_val.pos);
  valid_root = true;
  for (i = 0; i < NDIM; i++) {
    // initial root coords at first step are (0, 0, 0) and updated at each
    // iteration
    _xor[i] = xp[i] ^ Local.Root_Coords[i];
  }
  for (i = IMAX >> 1; i > root_val.level; i >>= 1) {
    for (j = 0; j < NDIM; j++) {
      if (_xor[j] & i) {
        valid_root = false;
        break;
      }
    }
    if (!valid_root) {
      break;
    }
  }
  if (!valid_root) {
    if (root != static_cast<cellptr>(G_root.get())) {
      root_level = Level((*root));
      for (j = i; j > root_level; j >>= 1) {
        // traverse up to the root level
        auto const parent = Parent(*root);
        root              = parent.get();
      }
      valid_root = true;
      root_val   = *root;
      for (i = IMAX >> 1; i > root_val.level; i >>= 1) {
        for (j = 0; j < NDIM; j++) {
          if (_xor[j] & i) {
            valid_root = false;
            break;
          }
        }
        if (!valid_root) {
          using const_ptr = decltype(bodytab)::const_iterator::pointer;
          printf(
              "P%ld body %ld\n",
              ProcessId,
              p - static_cast<const_ptr>(bodytab.begin()));
          root = G_root.get();
        }
      }
    }
  }
  root            = G_root.get();
  auto mynode     = static_cast<nodeptr>(root);
  auto mynode_val = static_cast<cell>(*NODE_AS_CELL(mynode));
  // determine which child holds the integerized coords of p
  kidIndex = subindex(xp, mynode_val.level);

  // get pointer to child
  auto *qptr = &(mynode_val.subp[kidIndex]);

  // go 1 level down and traverse the child
  l    = mynode_val.level >> 1;
  flag = true;

  /*
      Procedure Quad_Tree_Insert(j, n) … Try to insert particle j at node n in
     Quad_Tree
      if n an internal node              … n has 4 children
          CASE 1
          determine which child c of node n contains particle j
          Quad_Tree_Insert(j, c)
     else if n contains 1 particle   …  n is a leaf
          CASE 2
          add n’s 4 children to the Quad_Tree
          move the particle already in n into the child containing it
          let c be the child of n containing j
          Quad_Tree_Insert(j, c)
      else                                         …  n empty
          CASE 3
          store particle j in node n
      end
  */

  while (flag) { /* loop descending tree     */
    if (l == 0) {
      error("not enough levels in tree\n");
    }

    if (*qptr == nodeptr{}) {
      /* lock the parent cell */
      // CASE 3: Empty Cell!!
      //{pthread_mutex_lock(&CellLock->CL[((cellptr)mynode)->seqnum %
      // MAXLOCK]);}
      auto const sn = mynode_val.seqnum;
      CellLock.at(mynode_val.seqnum % MAXLOCK).lock();
      // reloading mynode_val since it may have changed after acquiring the
      // lock
      mynode_val = static_cast<cell>(*NODE_AS_CELL(mynode));
      ASSERT(mynode_val.seqnum == sn);
      if (*qptr == nodeptr{}) {
        // insert particle into the cell, treat the cell as a leaf
        le = InitLeaf(static_cast<cellptr>(mynode), ProcessId);
        ASSERT(le.local());
        auto &le_val = *(le.local());

        body_val.parent = le;
        body_val.level  = l;
        // Subindex of p in leaf
        body_val.child_num = le_val.num_bodies;
        // Subindex of leaf in cell
        le_val.child_num = kidIndex;
        // Assign body to leaf index
        le_val.bodyp[le_val.num_bodies++] = p;
        // Update current quad tree ptr to leaf
        flag = false;
        *p   = body_val;
        // This changes mynode_val
        *qptr = le;
        // write mynode_val back
        *NODE_AS_CELL(mynode) = mynode_val;
      }
      CellLock.at(mynode_val.seqnum % MAXLOCK).unlock();
      /* unlock the parent cell */
    }
    if (flag && (*qptr != nodeptr{}) && Type(*(*qptr)) == LEAF) {
      /*   reached a "leaf"?      */
      //{pthread_mutex_lock(&CellLock->CL[((cellptr)mynode)->seqnum %
      // MAXLOCK]);};
      auto const sn = mynode_val.seqnum;
      CellLock.at(mynode_val.seqnum % MAXLOCK).lock();
      /* lock the parent cell */
      // reloading mynode_val since it may have changed after acquiring the
      // lock
      mynode_val = static_cast<cell>(*NODE_AS_CELL(mynode));
      ASSERT(mynode_val.seqnum == sn);
      if (Type(*(*qptr)) == LEAF) { /* still a "leaf"?      */
        le          = *qptr;
        auto le_val = static_cast<leaf>(*le);
        if (le_val.num_bodies == MAX_BODIES_PER_LEAF) {
          // CASE 2: Subdivide the Tree
          // mynode_val.subp[kidIndex]) = SubdivideLeaf(le, mynode, l,
          // ProcessId)
          *qptr                 = SubdivideLeaf(le, mynode, l, ProcessId);
          *NODE_AS_CELL(mynode) = mynode_val;
        }
        else {
          // CASE 1: Leaf has still some free space, so add the particle

          body_val.parent                   = le;
          body_val.level                    = l;
          body_val.child_num                = le_val.num_bodies;
          le_val.bodyp[le_val.num_bodies++] = p;
          flag                              = false;
          // write it back
          *p  = body_val;
          *le = le_val;
        }
      }
      CellLock.at(mynode_val.seqnum % MAXLOCK).unlock();
      //{pthread_mutex_unlock(&CellLock->CL[((cellptr)mynode)->seqnum %
      // MAXLOCK]);};
      /* unlock the node           */
    }
    if (flag) {
      // We are not done yet (flag == true), so move down one level
      mynode = *qptr;
      // calculate and move down 1 level into that cell
      mynode_val = static_cast<cell>(*NODE_AS_CELL(mynode));
      kidIndex   = subindex(xp, l);
      qptr       = &(mynode_val.subp[kidIndex]); /* move down one level  */
      l          = l >> 1;                       /* and test next bit    */
    }
  }
  // update current local root
  SETV(Local.Root_Coords, xp);
  return Parent(*(*qptr));
}

/* * INTCOORD: compute integerized coordinates.  * Returns: true
unless rp was out of bounds.  */

/* integerized coordinate vector [0,IMAX) */
/* real coordinate vector (system coords) */
bool intcoord(long xp[NDIM], vector const rp)
{
  long   k;
  bool   inb;
  double xsc;

  inb = true;
  ASSERT(NDIM == 3);
  for (k = 0; k < NDIM; k++) {
    xsc = (rp[k] - g_rmin[k]) / g_rsize;
    if (0.0 <= xsc && xsc < 1.0) {
      xp[k] = floor(IMAX * xsc);
    }
    else {
      inb = false;
    }
  }
  return (inb);
}

/*
 * SUBINDEX: determine which subcell to select.
 */

/* integerized coordinates of particle */
/* current level of tree */
long subindex(long x[NDIM], long l)
{
  long i, k;
  bool yes;

  i   = 0;
  yes = false;
  if (x[0] & l) {
    i += NSUB >> 1;
    yes = true;
  }
  for (k = 1; k < NDIM; k++) {
    if (((x[k] & l) && !yes) || (!(x[k] & l) && yes)) {
      i += NSUB >> (k + 1);
      yes = true;
    }
    else
      yes = false;
  }

  return (i);
}

/*
 * HACKCOFM: descend tree finding center-of-mass coordinates.
 */

void hackcofm(long ProcessId)
{
  vector tmpv;

  /* get a cell using get*sub.  Cells are got in reverse of the order in */
  /* the cell array; i.e. reverse of the order in which they were created */
  /* this way, we look at child cells before parents       */

  for (auto ll = std::next(std::begin(Local.myleaftab), Local.mynleaf - 1);
       ll >= std::begin(Local.myleaftab);
       --ll) {
    auto l = *ll;
    ASSERT(l.local());
    auto &l_val = *l.local();
    l_val.mass  = 0.0;
    l_val.cost  = 0;
    CLRV(l_val.pos);
    for (size_t i = 0; i < static_cast<size_t>(l_val.num_bodies); i++) {
      auto const p_val = static_cast<body>(*(l_val.bodyp[i]));
      l_val.mass += p_val.mass;
      l_val.cost += p_val.cost;
      MULVS(tmpv, p_val.pos, p_val.mass);
      ADDV(l_val.pos, l_val.pos, tmpv);
    }
    DIVVS(l_val.pos, l_val.pos, l_val.mass);
#ifdef QUADPOLE
    CLRM(l_val.quad);
    for (size_t i = 0; i < l_val.num_bodies; i++) {
      p = l_val.bodyp[i];
      SUBV(dr, Pos(p), Pos(l));
      OUTVP(drdr, dr, dr);
      DOTVP(drsq, dr, dr);
      SETMI(Idrsq);
      MULMS(Idrsq, Idrsq, drsq);
      MULMS(tmpm, drdr, 3.0);
      SUBM(tmpm, tmpm, Idrsq);
      MULMS(tmpm, tmpm, Mass(p));
      ADDM(Quad(l), Quad(l), tmpm);
    }
#endif
    AtomicDone(*l) = true;
  }

  for (auto cc = std::next(std::begin(Local.mycelltab), Local.myncell - 1);
       cc >= std::begin(Local.mycelltab);
       --cc) {
    auto q = *cc;
    ASSERT(q.local());

    auto &q_val = *q.local();
    q_val.mass  = 0.0;
    q_val.cost  = 0;
    CLRV(q_val.pos);
    for (size_t i = 0; i < NSUB; i++) {
      auto r = q_val.subp[i];
      if (r != nodeptr{}) {
        while (!(static_cast<long>(AtomicDone(*r)))) {
          /* wait */
        }
        cell const r_val = *(NODE_AS_CELL(r));
        q_val.mass += r_val.mass;
        q_val.cost += r_val.cost;
        MULVS(tmpv, r_val.pos, r_val.mass);
        ADDV(q_val.pos, q_val.pos, tmpv);
        AtomicDone(*r) = false;
      }
    }
    DIVVS(q_val.pos, q_val.pos, q_val.mass);
#ifdef QUADPOLE
    // TODO rko: port to dash
    CLRM(Quad(q));
    for (i = 0; i < NSUB; i++) {
      r = Subp(q)[i];
      if (r != NULL) {
        SUBV(dr, Pos(r), Pos(q));
        OUTVP(drdr, dr, dr);
        DOTVP(drsq, dr, dr);
        SETMI(Idrsq);
        MULMS(Idrsq, Idrsq, drsq);
        MULMS(tmpm, drdr, 3.0);
        SUBM(tmpm, tmpm, Idrsq);
        MULMS(tmpm, tmpm, Mass(r));
        ADDM(tmpm, tmpm, Quad(r));
        ADDM(Quad(q), Quad(q), tmpm);
      }
    }
#endif
    AtomicDone(*q) = true;
  }
}

cellptr SubdivideLeaf(leafptr le, cellptr parent, long l, long ProcessId)
{
  long    i, index;
  long    xp[NDIM];
  bodyptr bodies[MAX_BODIES_PER_LEAF];
  long    num_bodies;

  /* first copy leaf's bodies to temp array, so we can reuse the leaf */
  auto le_val = static_cast<leaf>(*le);
  num_bodies  = le_val.num_bodies;
  /*
  for (i = 0; i < num_bodies; i++) {
    bodies[i]    = [i];
    Bodyp(le)[i] = NULL;
  }
  */
  std::copy(le_val.bodyp, le_val.bodyp + num_bodies, bodies);
  std::fill(le_val.bodyp, le_val.bodyp + num_bodies, bodyptr{nullptr});

  le_val.num_bodies = 0;
  /* create the parent cell for this subtree */
  auto c = InitCell(parent, ProcessId);
  ASSERT(c.local());
  auto &c_val = *(c.local());

  c_val.child_num = le_val.child_num;
  /* do first particle separately, so we can reuse le */
  auto p     = bodies[0];
  auto p_val = static_cast<body>(*p);
  intcoord(xp, p_val.pos);
  index             = subindex(xp, l);
  c_val.subp[index] = le;
  le_val.child_num  = index;
  le_val.parent     = c;
  le_val.level      = l >> 1;
  /* set stuff for body */
  p_val.parent    = le;
  p_val.child_num = le_val.num_bodies;
  p_val.level     = l >> 1;
  /* insert the body */
  le_val.bodyp[le_val.num_bodies++] = p;
  // write back body_val to guarantee consistency
  *p = p_val;
  // write back le_val to guarantee consistency
  *le = le_val;
  /* now handle the rest */
  for (i = 1; i < num_bodies; i++) {
    auto p = bodies[i];
    p_val  = static_cast<body>(*p);
    intcoord(xp, p_val.pos);
    index = subindex(xp, l);
    if (c_val.subp[index] == nodeptr{}) {
      le                = InitLeaf(c, ProcessId);
      le_val            = *le;
      le_val.child_num  = index;
      c_val.subp[index] = le;
    }
    else {
      le     = c_val.subp[index];
      le_val = *le;
    }
    p_val.parent                      = le;
    p_val.child_num                   = le_val.num_bodies;
    p_val.level                       = l >> 1;
    le_val.bodyp[le_val.num_bodies++] = p;
    // write it back to guarantee consistency
    *p  = p_val;
    *le = le_val;
  }

  return c;
}

/*
 * MAKECELL: allocation routine for cells.
 */

cellptr makecell(long ProcessId)
{
  if (Local.mynumcell == maxmycell) {
    error(
        "makecell: Proc %ld needs more than %ld cells; increase fcells\n",
        ProcessId,
        maxmycell);
  }
  auto const Mycell       = Local.mynumcell++;
  auto  c_gpos = celltab.pattern().global(Mycell);
  auto  c_gptr = static_cast<cellptr>(celltab.begin() + c_gpos);
  auto &c      = celltab.local[Mycell];
  c.seqnum     = ProcessId * maxmycell + Mycell;
  c.type       = CELL;
  c.done       = false;
  c.mass       = 0.0;
  c.cost       = 0;

  std::fill(std::begin(c.pos), std::end(c.pos), real{0});
  /*
  for (i = 0; i < NSUB; i++) {
    Subp(c)[i] = NULL;
  }
  */
  std::fill(c.subp, c.subp + NSUB, cellptr{nullptr});

  Local.mycelltab[Local.myncell++] = c_gptr;
  return c_gptr;
}

/*
 * MAKELEAF: allocation routine for leaves.
 */

leafptr makeleaf(long ProcessId)
{
  if (Local.mynumleaf == maxmyleaf) {
    error(
        "makeleaf: Proc %ld needs more than %ld leaves; increase fleaves\n",
        ProcessId,
        maxmyleaf);
  }

  auto const Myleaf = Local.mynumleaf++;
  // Resolve the global index of a specific unit's local index 'Myleaf'
  auto const gpos = leaftab.pattern().global(Myleaf);
  // Get global pointer to local element
  auto le_ptr = static_cast<leafptr>(leaftab.begin() + gpos);
  ASSERT(le_ptr.local());
  auto &le_val      = leaftab.local[Myleaf];
  le_val.seqnum     = ProcessId * maxmyleaf + Myleaf;
  le_val.type       = LEAF;
  le_val.done       = false;
  le_val.mass       = 0.0;
  le_val.num_bodies = 0;
  le_val.cost       = 0;

  std::fill(le_val.bodyp, le_val.bodyp + MAX_BODIES_PER_LEAF, bodyptr{});
  std::fill(std::begin(le_val.pos), std::end(le_val.pos), real{0});

  Local.myleaftab[Local.mynleaf++] = le_ptr;
  return (le_ptr);
}
