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

#include "code.h"
#include "load.h"

extern pthread_t PThreadTable[];

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
  Local.Current_Root = G_root.get().get();

  for (bodyptr *pp = Local.mybodytab; pp < Local.mybodytab + Local.mynbody;
       ++pp) {
    bodyptr p        = *pp;
    auto const pmass = Mass(*p).get();
    // insert all bodies into the tree
    if (pmass != 0.0) {
      Local.Current_Root =
          loadtree(p, static_cast<cellptr>(Local.Current_Root), ProcessId);
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
  hackcofm(ProcessId);
  dash::barrier();
}

cellptr InitCell(cellptr parent, long ProcessId)
{
  auto c          = makecell(ProcessId);
  auto c_val      = static_cast<cell>(*c);
  c_val.processor = ProcessId;
  c_val.next      = NULL;
  c_val.prev      = NULL;
  if (parent == cellptr{nullptr})
    c_val.level = IMAX >> 1;
  else
    c_val.level = Level((*parent)) >> 1;

  c_val.parent    = parent;
  c_val.child_num = 0;
  *c              = c_val;
  return c;
}

leafptr InitLeaf(cellptr parent, long ProcessId)
{
  auto l     = makeleaf(ProcessId);
  auto l_ref = *l;
  // dereference it, local copy
  auto lv      = static_cast<leaf>(l_ref);
  lv.processor = ProcessId;
  lv.next      = NULL;
  lv.prev      = NULL;
  if (parent == cellptr{nullptr})
    lv.level = IMAX >> 1;
  else
    lv.level = Level((*parent)) >> 1;
  lv.parent  = parent;

  lv.child_num = 0;

  // write it back
  l_ref = lv;
  return (l);
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
        if (!c_val.subp[k]) {
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
        if (c_val.subp[k]) {
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
  std::cout << ss.str() << std::endl;
}

/*
 * LOADTREE: descend tree and insert particle.
 */

nodeptr loadtree(bodyptr p, cellptr root, long ProcessId)
{
  long l, xp[NDIM], _xor[NDIM], flag;
  long i, j, root_level;
  bool valid_root;
  long kidIndex;
  nodeptr *qptr, mynode;
  cellptr mycell;
  leafptr le;

  auto body_val = static_cast<body>(*p);
  auto root_val = static_cast<cell>(*root);
  // get integerized coords of real position in space
  intcoord(xp, body_val.pos);
  valid_root = TRUE;
  for (i = 0; i < NDIM; i++) {
    // initial root coords at first step are (0, 0, 0) and update at each
    // iteration
    _xor[i] = xp[i] ^ Local.Root_Coords[i];
  }
  for (i = IMAX >> 1; i > root_val.level; i >>= 1) {
    for (j = 0; j < NDIM; j++) {
      if (_xor[j] & i) {
        valid_root = FALSE;
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
        dash::GlobRef<nodeptr> parent = Parent((*root));
        cellptr parent_cell           = static_cast<nodeptr>(parent);
        root                          = parent_cell;
      }
      valid_root = TRUE;
      root_val   = *root;
      for (i = IMAX >> 1; i > root_val.level; i >>= 1) {
        for (j = 0; j < NDIM; j++) {
          if (_xor[j] & i) {
            valid_root = FALSE;
            break;
          }
        }
        if (!valid_root) {
          using const_ptr = decltype(bodytab)::const_iterator::pointer;
          printf("P%ld body %ld\n", ProcessId,
                 p - static_cast<const_ptr>(bodytab.begin()));
          root = G_root.get();
        }
      }
    }
  }
  root            = G_root.get();
  mynode          = root;
  auto mynode_val = static_cast<cell>(*NODE_AS_CELL(mynode));
  // determine which child holds the integerized coords of p
  kidIndex = subindex(xp, mynode_val.level);

  // get pointer to child
  qptr = &(mynode_val.subp[kidIndex]);

  // go 1 level down and traverse the child
  l    = mynode_val.level >> 1;
  flag = TRUE;

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

    // verify if we have a null pointer
    constexpr auto const nullgptr = nodeptr{};
    if (*qptr == nullgptr) {
      /* lock the parent cell */
      // CASE 3: Empty Cell!!
      //{pthread_mutex_lock(&CellLock->CL[((cellptr)mynode)->seqnum %
      // MAXLOCK]);}
      CellLock.at(mynode_val.seqnum % MAXLOCK).lock();
      if (*qptr == nullgptr) {
        // insert particle into the cell, treat the cell as a leaf
        le = InitLeaf(static_cast<cellptr>(mynode), ProcessId);
        auto le_val = static_cast<leaf>(*le);

        body_val.parent = le;
        body_val.level  = l;
        // Subindex of p in leaf
        body_val.child_num = le_val.num_bodies;
        // Subindex of leaf in cell
        le_val.child_num = kidIndex;
        // Assign body to leaf index
        le_val.bodyp[le_val.num_bodies++] = p;
        // Update current quad tree ptr to leaf
        flag = FALSE;
        *p   = body_val;
        *le  = le_val;
        // mynode_val.subp[kidIndex]) = le
        *qptr = le;
        // cast as a cell
        *NODE_AS_CELL(mynode) = mynode_val;
#ifdef ENABLE_ASSERTIONS
        auto const tmp        = static_cast<cell>(*NODE_AS_CELL(mynode));
        ASSERT(tmp.subp[kidIndex] == le);
#endif
      }
      CellLock.at(mynode_val.seqnum % MAXLOCK).unlock();
      /* unlock the parent cell */
    }
    if (flag && (*qptr) && Type(*(*qptr)) == LEAF) {
      /*   reached a "leaf"?      */
      //{pthread_mutex_lock(&CellLock->CL[((cellptr)mynode)->seqnum %
      // MAXLOCK]);};
      CellLock.at(mynode_val.seqnum % MAXLOCK).lock();
      /* lock the parent cell */
      if (Type(*(*qptr)) == LEAF) { /* still a "leaf"?      */
        le          = *qptr;
        auto le_val = static_cast<leaf>(*le);
        if (le_val.num_bodies == MAX_BODIES_PER_LEAF) {
          // CASE 2: Subdivide the Tree
          // mynode_val.subp[kidIndex]) = SubdivideLeaf(le, mynode, l,
          // ProcessId)
          // TODO: is this really correct or does it not loose some updates to
          // mynode in SubdivideLeaf
          *qptr                 = SubdivideLeaf(le, mynode, l, ProcessId);
          *NODE_AS_CELL(mynode) = mynode_val;
        }
        else {
          // CASE 1: Leaf has still some free space, so add the particle

          body_val.parent                   = le;
          body_val.level                    = l;
          body_val.child_num                = le_val.num_bodies;
          le_val.bodyp[le_val.num_bodies++] = p;
          flag                              = FALSE;
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
  nodeptr tmp_ptr = *qptr;
  return Parent((*tmp_ptr));
}

/* * INTCOORD: compute integerized coordinates.  * Returns: TRUE
unless rp was out of bounds.  */

/* integerized coordinate vector [0,IMAX) */
/* real coordinate vector (system coords) */
bool intcoord(long xp[NDIM], vector const rp)
{
  long k;
  bool inb;
  double xsc;

  inb = TRUE;
  ASSERT(NDIM == 3);
  sh_vec rmin_val = rmin.get();
  vector rmin_tmp = {rmin_val.x, rmin_val.y, rmin_val.z};
  for (k = 0; k < NDIM; k++) {
    xsc = (rp[k] - rmin_tmp[k]) / rsize.get();
    if (0.0 <= xsc && xsc < 1.0) {
      xp[k] = floor(IMAX * xsc);
    }
    else {
      inb = FALSE;
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
  long yes;

  i   = 0;
  yes = FALSE;
  if (x[0] & l) {
    i += NSUB >> 1;
    yes = TRUE;
  }
  for (k = 1; k < NDIM; k++) {
    if (((x[k] & l) && !yes) || (!(x[k] & l) && yes)) {
      i += NSUB >> (k + 1);
      yes = TRUE;
    }
    else
      yes = FALSE;
  }

  return (i);
}

/*
 * HACKCOFM: descend tree finding center-of-mass coordinates.
 */

void hackcofm(long ProcessId)
{
  long i;
  nodeptr r;
  leafptr l;
  leafptr *ll;
  bodyptr p;
  cellptr q;
  cellptr *cc;
  vector tmpv;

  /* get a cell using get*sub.  Cells are got in reverse of the order in */
  /* the cell array; i.e. reverse of the order in which they were created */
  /* this way, we look at child cells before parents       */

  dash::GlobRef<cellptr> cellptr_ref = G_root.get();
  cellptr cellptr_val = cellptr_ref.get();

  for (ll = Local.myleaftab + Local.mynleaf - 1; ll >= Local.myleaftab;
       ll--) {
    l          = *ll;
    auto l_val = static_cast<leaf>(*l);
    l_val.mass = 0.0;
    l_val.cost = 0;
    CLRV(l_val.pos);
    for (i = 0; i < l_val.num_bodies; i++) {
      p                = l_val.bodyp[i];
      auto const p_val = static_cast<body>(*p);
      l_val.mass += p_val.mass;
      l_val.cost += p_val.cost;
      MULVS(tmpv, p_val.pos, p_val.mass);
      ADDV(l_val.pos, l_val.pos, tmpv);
    }
    DIVVS(l_val.pos, l_val.pos, l_val.mass);
#ifdef QUADPOLE
    // TODO: port to dash
    CLRM(l_val.quad);
    for (i = 0; i < l_val.num_bodies; i++) {
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
    // write the modified value back
    *l = l_val;
    // write atomically true
    AtomicDone(*l) = TRUE;
  }

  for (cc = Local.mycelltab + Local.myncell - 1; cc >= Local.mycelltab;
       cc--) {
    q          = *cc;

    if (cellptr_val == q) {
//      std::cout << "I am here" << std::endl;
    }
    auto q_val = static_cast<cell>(*q);
    q_val.mass = 0.0;
    q_val.cost = 0;
    CLRV(q_val.pos);
    for (i = 0; i < NSUB; i++) {
      r = q_val.subp[i];
      if (r) {
        while (!(static_cast<long>(AtomicDone(*r)))) {
          /* wait */
        }
        cell const r_val = *(NODE_AS_CELL(r));
        //std::cout << "rval:\n" << r_val << std::endl;
        q_val.mass += r_val.mass;
        q_val.cost += r_val.cost;
        MULVS(tmpv, r_val.pos, r_val.mass);
        //std:: cout << "[" << tmpv[0] << ", " << tmpv[1] << ", " << tmpv[2] << "]" << std::endl;
        ADDV(q_val.pos, q_val.pos, tmpv);
    if (cellptr_val == q) {
      //std::cout << "Value of q_val is::\n" << q_val << std::endl;
    }
        AtomicDone(*r) = FALSE;
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
    *q             = q_val;
    if (cellptr_val == q) {
      //std::cout << "Value of qval before write back is:\n" << q_val << std::endl;
    }
    AtomicDone(*q) = TRUE;
  }
}

cellptr SubdivideLeaf(leafptr le, cellptr parent, long l, long ProcessId)
{
  long i, index;
  long xp[NDIM];
  bodyptr bodies[MAX_BODIES_PER_LEAF];
  long num_bodies;

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
  auto c     = InitCell(parent, ProcessId);
  auto c_val = static_cast<cell>(*c);

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
  /* now handle the rest */
  // write back body_val to guarantee consistency
  *p = p_val;
  // write back le_val to guarantee consistency
  *le = le_val;
  for (i = 1; i < num_bodies; i++) {
    p     = bodies[i];
    p_val = static_cast<body>(*p);
    intcoord(xp, p_val.pos);
    index = subindex(xp, l);
    if (!c_val.subp[index]) {
      le                = InitLeaf(c, ProcessId);
      le_val            = *le;
      le_val.child_num  = index;
      c_val.subp[index] = le;
    }
    else {
      le = c_val.subp[index];
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
  *c = c_val;
  return c;
}

/*
 * MAKECELL: allocation routine for cells.
 */

cellptr makecell(long ProcessId)
{
  long Mycell;

  if (Local.mynumcell == maxmycell) {
    error("makecell: Proc %ld needs more than %ld cells; increase fcells\n",
          ProcessId, maxmycell);
  }
  Mycell       = Local.mynumcell++;
  auto c_gpos  = celltab.pattern().global(Mycell);
  auto c       = static_cast<cellptr>(celltab.begin() + c_gpos);
  auto c_val   = static_cast<cell>(*c);
  c_val.seqnum = ProcessId * maxmycell + Mycell;
  c_val.type   = CELL;
  // TODO: verify if we need atomic here
  c_val.done = FALSE;
  c_val.mass = 0.0;
  /*
  for (i = 0; i < NSUB; i++) {
    Subp(c)[i] = NULL;
  }
  */
  std::fill(c_val.subp, c_val.subp + NSUB, cellptr{nullptr});

  *c = c_val;

  Local.mycelltab[Local.myncell++] = c;
  return c;
}

/*
 * MAKELEAF: allocation routine for leaves.
 */

leafptr makeleaf(long ProcessId)
{
  long Myleaf;

  if (Local.mynumleaf == maxmyleaf) {
    error("makeleaf: Proc %ld needs more than %ld leaves; increase fleaves\n",
          ProcessId, maxmyleaf);
  }

  Myleaf          = Local.mynumleaf++;
  //Resolve the global index of a specific unit's local index 'Myleaf'
  auto const gpos = leaftab.pattern().global(Myleaf);
  auto le         = static_cast<leafptr>(leaftab.begin() + gpos);
  auto le_val     = static_cast<leaf>(*le);
  le_val.seqnum   = ProcessId * maxmyleaf + Myleaf;
  le_val.type     = LEAF;
  // TODO: verify if we need atomic here
  le_val.done       = FALSE;
  le_val.mass       = 0.0;
  le_val.num_bodies = 0;

  std::fill(le_val.bodyp, le_val.bodyp + MAX_BODIES_PER_LEAF, bodyptr{});
  // write the value back to global memory
  *le                              = le_val;
  Local.myleaftab[Local.mynleaf++] = le;
  return (le);
}
