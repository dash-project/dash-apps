#ifndef COMM_H_INCLUDED
#define COMM_H_INCLUDED

#include <array>
#include "mpi.h"
#include "domain.h"
#include "lulesh-mock.h"

#define MAX_FIELDS_PER_MPI_COMM  6

class Domain;
typedef Real_t& (Domain::*Domain_member)(Index_t);


//
// encapsulate all communication functionality
// 
class Comm
{
 public:
  // what is being transferred
  enum Transfer 
  {
    X0=0, X1, Y0, Y1,        // edges
    X0Y0, X0Y1, X1Y0, X1Y1,  // corners
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
  std::array<int, 2> m_nProcs;
    
  Transfer    m_trans[3][3];
  int         m_neigh[3][3];
  MPI_Request m_reqs[3][3];

  std::vector<Real_t> m_sendBuf;
  std::vector<Real_t> m_recvBuf;
  
  std::array<size_t, 8> m_bufOffset;
  
 public:
  Comm(Domain& dom, std::array<int,2> nProcs);

  Index_t numProcs() const { return m_numProcs; }
  Index_t myRank() const { return m_myRank; }

  // covert (col, row) to rank
  Index_t toRank(Index_t col, Index_t row) {
    if( col<0 || col>=px() )
      return -1;
    if( row<0 || row>=py() )
      return -1;

   // WORKAROUND because DASH currently does not 
   // allow COL_MAJOR teamspec: 
   // row/col and px() py() are switched 
    return col*py()+row;
  }
  
  // number of processes in x and y direction
  Index_t px() const { return m_nProcs[0]; }
  Index_t py() const { return m_nProcs[1]; }

   // WORKAROUND because DASH currently does not 
   // allow COL_MAJOR teamspec: 
   // row/col and px() py() are switched 

  Index_t col() const { return myRank()/py(); }
  Index_t row() const { return myRank()%py(); }
  
  Real_t* sendBuf(Transfer trans) {
    return m_sendBuf.data()+m_bufOffset[trans];
  }
  
  Real_t* recvBuf(Transfer trans) {
    return m_recvBuf.data()+m_bufOffset[trans];
  }

  size_t MsgSize(std::array<int, 2> dim,
		 Transfer trans,
		 std::vector<Domain_member>& fields);
  
  void Send(std::array<int, 2> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);

  void Recv(std::array<int, 2> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);
  
  void Sync(std::array<int, 2> dim,
	    std::vector<Domain_member>& fields,
	    Scope scope,
	    Direction dir,
	    Action action);

 private:
  
  //
  // packing and unpacking of data elements between the domain data
  // structure and the provided buffer
  //
  void PackUnpack(std::array<int, 2> dim,
		  std::vector<Domain_member>& fields,
		  Transfer trans,
		  Real_t *buffer,
		  Pack pack,
		  Action action);
    
};

#endif /* COMM_H_INCLUDED */

