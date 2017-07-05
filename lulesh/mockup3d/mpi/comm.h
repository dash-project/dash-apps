#ifndef COMM_H_INCLUDED
#define COMM_H_INCLUDED

#include <array>
#include "mpi.h"
#include "domain.h"
#include "lulesh-mock.h"

#define MAX_FIELDS_PER_MPI_COMM  6

class Domain;
typedef Real_t& (Domain::*Domain_member)(Index_t);

#define IS_PLANE(trans_)			\
  (trans_==X0 || trans_==X1 || trans_==Y0 ||	\
   trans_==Y1 || trans_==Z0 || trans_==Z1)

#define IS_EDGE(trans_)							\
  (trans_==X0Y0 || trans_==X0Y1 || trans_==X1Y0 || trans_==X1Y1 ||	\
   trans_==X0Z0 || trans_==X0Z1 || trans_==X1Z0 || trans_==X1Z1 ||	\
   trans_==Y0Z0 || trans_==Y0Z1 || trans_==Y1Z0 || trans_==Y1Z1)

#define IS_CORNER(trans_)						\
  (trans_==X0Y0Z0 || trans_==X0Y1Z0 || trans_==X1Y0Z0 || trans_==X1Y1Z0 || \
   trans_==X0Y0Z1 || trans_==X0Y1Z1 || trans_==X1Y0Z1 || trans_==X1Y1Z1)

//
// encapsulate all communication functionality
// 
class Comm
{
 public:
  // what is being transferred
  enum Transfer 
  {
    X0=0, X1, Y0, Y1, Z0, Z1,  // 6 planes
    X0Y0, X0Y1, X1Y0, X1Y1,    // 12 edges
    X0Z0, X0Z1, X1Z0, X1Z1,
    Y0Z0, Y0Z1, Y1Z0, Y1Z1,
    X0Y0Z0, X0Y0Z1, X0Y1Z0, X0Y1Z1,  // 8 corners
    X1Y0Z0, X1Y0Z1, X1Y1Z0, X1Y1Z1,
    NONE
  };
  
  // in which direction the communication happens
  enum class Direction { One, Both };

  // what to do with incoming values
  enum class Action { Add, Replace };

  // what to communicate
  // in 3D: All=corners+edges+planes, Plane=planes
  // in 2D: All=corners+edges, Plane=edges
  enum class Scope { All, Plane };

  // toggle between packing and unpacking
  enum class Pack { Yes, No };
  
 private:
  Domain& m_dom;

  Index_t m_numProcs;
  Index_t m_myRank;

  // number of processes in each dim.
  std::array<int, 3> m_nProcs;
    
  Transfer    m_trans[3][3][3];
  int         m_neigh[3][3][3];
  MPI_Request m_reqs[3][3][3];

  std::vector<Real_t> m_sendBuf;
  std::vector<Real_t> m_recvBuf;
  
  std::array<size_t, 26> m_bufOffset;
  
 public:
  Comm(Domain& dom, std::array<int,3> nProcs);

  Index_t numProcs() const { return m_numProcs; }
  Index_t myRank() const { return m_myRank; }

  // covert (col, row, plane) to rank
  Index_t toRank(Index_t col, Index_t row, Index_t plane) {
    if( col<0 || col>=px() )
      return -1;
    if( row<0 || row>=py() )
      return -1;
    if( plane<0 || plane>=pz() )
      return -1;

    return plane*px()*py()+row*px()+col;
  }
  
  // number of processes in x and y direction
  Index_t px() const { return m_nProcs[0]; }
  Index_t py() const { return m_nProcs[1]; }
  Index_t pz() const { return m_nProcs[2]; }
  
  Index_t plane() const { return myRank()/(px()*py()); }
  Index_t row()   const { return (myRank()-plane()*px()*py())/px(); }
  Index_t col()   const { return (myRank()-plane()*px()*py())%px(); }

  size_t offset(Transfer trans) {
    return m_bufOffset[trans];
  }
  
  Real_t* sendBuf(Transfer trans) {
    return m_sendBuf.data()+offset(trans);
  }
  
  Real_t* recvBuf(Transfer trans) {
    return m_recvBuf.data()+offset(trans);
  }

  size_t MsgSize(std::array<int,3> dim,
		 Transfer trans,
		 std::vector<Domain_member>& fields);
  
  void Send(std::array<int,3> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);

  void Recv(std::array<int,3> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);
  
  void Sync(std::array<int,3> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);

 private:
  
  //
  // packing and unpacking of data elements between the domain data
  // structure and the provided buffer
  //
  void PackUnpack(std::array<int,3> dim,
		  std::vector<Domain_member>& fields,
		  Transfer trans,
		  Real_t *buffer,
		  Pack pack,
		  Action action);
    
};

#endif /* COMM_H_INCLUDED */

