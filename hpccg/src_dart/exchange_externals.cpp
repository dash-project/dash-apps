
//@HEADER
// ************************************************************************
// 
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifdef USING_MPI  // Compile this routine only if running in parallel
#include <iostream>
using std::cerr;
using std::endl;
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <libdash.h>
#include "exchange_externals.hpp"
#undef DEBUG
void exchange_externals(HPC_Sparse_Matrix * A, const double *x)
{
  int i, j, k;
  int num_external = 0;

  // Extract Matrix pieces

  int local_nrow = A->local_nrow;
  int num_neighbors = A->num_send_neighbors;
  int * recv_length = A->recv_length;
  int * send_length = A->send_length;
  int * neighbors = A->neighbors;
  double * send_buffer = A->send_buffer;
  int total_to_be_sent = A->total_to_be_sent;
  int * elements_to_send = A->elements_to_send;

  dart_team_unit_t myid;
  dart_team_myid(DART_TEAM_ALL, &myid);

  //
  //  first post receives, these are immediate receives
  //  Do not wait for result to come, will do that at the
  //  wait call below.
  //

//  int MPI_MY_TAG = 99;
//
//  MPI_Request * request = new MPI_Request[num_neighbors];

  //
  // Externals are at end of locals
  //
  double *x_external = (double *) x + local_nrow;

//  // Post receives first
//  for (i = 0; i < num_neighbors; i++)
//    {
//      int n_recv = recv_length[i];
//      MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
//		MPI_COMM_WORLD, request+i);
//      x_external += n_recv;
//    }


  //
  // Fill up send buffer
  //

  for (i=0; i<total_to_be_sent; i++) send_buffer[i] = x[elements_to_send[i]];

  //
  // Send to each neighbor
  //

//  std::cout << "["<< myid.id << "] " << "Communicating with " << num_neighbors << " neighbors" << std::endl;

  for(int i=0; i<num_neighbors; ++i) {
    int n_send = send_length[i];
//    std::cout << "["<< myid.id << "] " << "Putting " << n_send << " elements to gptr " << A->gptr[i] << std::endl;
    dart_put(A->gptr[i], send_buffer, n_send, DART_TYPE_DOUBLE);
    send_buffer += n_send;
  }

  dart_flush_all(A->data_win.gptr);

  int val = 1;
  // signal completion
  dart_gptr_t gptr = A->signal_win.gptr;
  for(int i=0; i<num_neighbors; ++i) {
    dart_gptr_setunit(&gptr, DART_TEAM_UNIT_ID(neighbors[i]));
//    std::cout << "["<< myid.id << "] " << "Sending completion flag to " << neighbors[i] << std::endl;
    dart_accumulate(gptr, &val, 1, DART_TYPE_INT, DART_OP_SUM);
  }

  // ... and wait for completion of neighbors
  int flag = 0;
  dart_gptr_setunit(&gptr, myid);
  do {
    const int zero = 0;

//    dart_fetch_and_op(gptr, &zero, &val, DART_TYPE_INT, DART_OP_NO_OP);
//    std::cout << "["<< myid.id << "] value: " << val << std::endl;
//    std::cout << "["<< myid.id << "] gptr: " << gptr << std::endl;
    dart_compare_and_swap(gptr, &zero, &num_neighbors, &flag, DART_TYPE_INT);
    dart_flush_all(gptr);
//    std::cout << "["<< myid.id << "] flag: " << flag << std::endl;
//    std::cout << "["<< myid.id << "] baseptr (" << A->signal_win.baseptr
//              << "): " << *static_cast<int*>(A->signal_win.baseptr)
//              << std::endl;
//    sleep(1);
  } while (flag != num_neighbors);
//  std::cout << "["<< myid.id << "]All puts seem to be completed!" << std::endl;


  // copy externals to end of vector
  double *begin = reinterpret_cast<double*>(A->data_win.baseptr);
  double *end   = begin + A->total_received_elements;
//  std::cout << "["<< myid.id << "] Copying " << A->total_received_elements << " from "  << begin << " to " << x_external << std::endl;
  std::copy(begin, end, x_external);

//  for (i = 0; i < num_neighbors; i++)
//    {
//      int n_send = send_length[i];
//      MPI_Send(send_buffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG,
//	       MPI_COMM_WORLD);
//      send_buffer += n_send;
//    }
//
//  //
//  // Complete the reads issued above
//  //
//
//  MPI_Status status;
//  for (i = 0; i < num_neighbors; i++)
//    {
//      if ( MPI_Wait(request+i, &status) )
//	{
//	  cerr << "MPI_Wait error\n"<<endl;
//	  exit(-1);
//	}
//    }
//
//  delete [] request;

  return;
}
#endif // USING_MPI
