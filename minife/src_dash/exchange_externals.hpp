#ifndef _exchange_externals_hpp_
#define _exchange_externals_hpp_

//@HEADER
// ************************************************************************
//
// MiniFE: Simple Finite Element Assembly and Solve
// Copyright (2006-2013) Sandia Corporation
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
//
// ************************************************************************
//@HEADER

#include <cstdlib>
#include <iostream>
#include <algorithm>

#ifdef HAVE_MPI
#include <mpi.h>
#include <dash/dart/if/dart.h>
#include <libdash.h>
#endif

#include <outstream.hpp>

#include <TypeTraits.hpp>

namespace miniFE {

template<typename MatrixType,
         typename VectorType>
void
finish_exchange_externals(MatrixType& A, VectorType& x);

template<typename MatrixType,
         typename VectorType>
void
exchange_externals(MatrixType& A,
                   VectorType& x)
{
#ifdef HAVE_MPI
#ifdef MINIFE_DEBUG
  std::ostream& os = outstream();
  os << "entering exchange_externals\n";
#endif

  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if (numprocs < 2) return;

  // issue async puts
  begin_exchange_externals(A, x);

  // use the finish for async here
  finish_exchange_externals(A, x);

#ifdef MINIFE_DEBUG
  os << "leaving exchange_externals"<<std::endl;
#endif

//endif HAVE_MPI
#endif
}

#ifdef HAVE_MPI
static std::vector<MPI_Request> exch_ext_requests;
#endif

template<typename MatrixType,
         typename VectorType>
void
begin_exchange_externals(MatrixType& A,
                         VectorType& x)
{
//  std::cout << "entering begin_exchange_externals\n";
#ifdef HAVE_MPI


  int numprocs = 1, myproc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  if (numprocs < 2) return;

  typedef typename MatrixType::ScalarType Scalar;
  typedef typename MatrixType::LocalOrdinalType LocalOrdinal;
  typedef typename MatrixType::GlobalOrdinalType GlobalOrdinal;

  // Extract Matrix pieces

  int num_neighbors = A.neighbors.size();
  const std::vector<LocalOrdinal>&  recv_length = A.recv_length;
  const std::vector<LocalOrdinal>&  send_length = A.send_length;
  const std::vector<int>&           neighbors   = A.neighbors;
  const std::vector<GlobalOrdinal>& elements_to_send = A.elements_to_send;

  std::vector<Scalar>& send_buffer = A.send_buffer;

  dart_datatype_t dart_dtype = TypeTraits<Scalar>::dart_type();

  //
  // Fill up send buffer
  //

  size_t total_to_be_sent = elements_to_send.size();
  for(size_t i=0; i<total_to_be_sent; ++i) {
    send_buffer[i] = x.coefs[elements_to_send[i]];
  }

  //
  // Send to each neighbor
  //

  Scalar* s_buffer = &send_buffer[0];

  for(int i=0; i<num_neighbors; ++i) {
    int n_send = send_length[i];
    int neighbor = i;
    int offset_at_neighbor = A.offset[i];
    // operator<< seems broken
//    auto out_it = A.data.row(i).begin() + A.offset[i];
//    std::cout << out_it << std::endl;
    std::cout << "copy_async to " << i << " at offset " << A.offset[i] << std::endl;
    // ends up in a DART_GPTR_NULL being passed to dart_put
    dash::copy_async(s_buffer, s_buffer + n_send + 1, A.data.row(i).begin() + offset_at_neighbor);
    s_buffer += n_send;
  }
#endif
}

template<typename MatrixType,
         typename VectorType>
void
finish_exchange_externals(MatrixType& A, VectorType& x)
{
#ifdef HAVE_MPI

  //
  // Complete the puts issued above
  //

  dart_team_unit_t myid;
  dart_team_myid(DART_TEAM_ALL, &myid);
  const int num_neighbors = A.neighbors.size();
  // complete put operations
  A.data.flush_all();

  const std::vector<int>& neighbors = A.neighbors;
  // signal completion
  for(int i=0; i<num_neighbors; ++i) {
    A.signal[neighbors[i]].add(1);
  }
  auto sig = A.signal[myid.id];
  while (sig.compare_exchange(num_neighbors, 0)) {};

  typedef typename MatrixType::ScalarType Scalar;
  // copy externals to vector
  std::vector<Scalar>& x_coefs = x.coefs;
  int local_nrow = A.rows.size();
  Scalar* x_external = &(x_coefs[local_nrow]);
  Scalar *begin = A.data.local.lbegin();
  Scalar *end   = begin + A.total_received_elements;
  std::copy(begin, end, x_external);
//endif HAVE_MPI
#endif
}

}//namespace miniFE

#endif

