
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
  // Extract Matrix pieces
  int local_nrow = A->local_nrow;
  int num_neighbors = A->num_send_neighbors;
  int * send_length = A->send_length;
  int * neighbors = A->neighbors;
  double * send_buffer = A->send_buffer;
  int total_to_be_sent = A->total_to_be_sent;
  int * elements_to_send = A->elements_to_send;

  //
  // Externals are at end of locals
  //
  double *x_external = (double *) x + local_nrow;

  //
  // Fill up send buffer
  //

  for (int i=0; i<total_to_be_sent; i++)
    send_buffer[i] = x[elements_to_send[i]];

  //
  // Initiate transfer to each neighbor
  //

  auto& narray = A->data;
  for(int i=0; i<num_neighbors; ++i) {
    int n_send = send_length[i];
    int neighbor = neighbors[i];
    int offset_at_neighbor = A->offset[i];
    dash::copy_async(send_buffer, send_buffer + n_send + 1,
      narray(neighbor).begin() + A->offset[i]);
    send_buffer += n_send;
  }

  // wait for the transfer to complete
  narray.flush_all();

  // signal completion to neighors
  for(int i=0; i<num_neighbors; ++i) {
    A->signal(neighbors[i]).post();
  }

  // ... and wait for completion of neighbors
  A->signal.wait();

  // copy externals to end of vector
  double *begin = A->data.lbegin();
  double *end   = begin + A->total_received_elements;
//  std::cout << "["<< myid.id << "] Copying " << A->total_received_elements << " from "  << begin << " to " << x_external << std::endl;
  std::copy(begin, end, x_external);


  return;
}
#endif // USING_MPI
