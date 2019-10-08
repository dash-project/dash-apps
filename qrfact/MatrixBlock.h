#ifndef CHOLESKY_MATRIXBLOCK_H
#define CHOLESKY_MATRIXBLOCK_H

#include <libdash.h>

//#define DEBUG

/**
 * Implementation of a matrix block. The block lies in a two-dimensional matrix
 * represented by a 4-dimensional DASH matrix where the first two dimensions
 * represent superblocks while the last two dimensions are the blocks inside the super block.
 * The MatrixBlock abstracts away the two superblock dimensions.
 */

template<typename MatrixT>
class MatrixBlock {
public:

  using value_t = typename MatrixT::value_type;

  MatrixBlock(
    MatrixT &matrix,
    int block_row_idx,
    int block_col_idx)
  : _matrix(&matrix) {
    auto& pattern       = matrix.pattern();
    int block_size   = pattern.extent(2); // blocks are in dim 2 and 3
    std::array<typename MatrixT::pattern_type::index_type, 4> gcoords = {block_row_idx, block_col_idx, 0, 0};
    this->_glob_idx = pattern.global_at(gcoords);
    typename MatrixT::pattern_type::local_index_t local_idx = pattern.local_index(gcoords);
#ifdef DEBUG_OUTPUT
    std::cout << "gcoords: {" << block_row_idx << "," << block_col_idx << "}, glob_idx: "<< _glob_idx << ", local_idx: {u:" << local_idx.unit.id << ", off: " << local_idx.index << "}" << std::endl;
#endif // DEBUG_OUTPUT
    this->_is_local = (local_idx.unit.id == dash::myid());
    this->_size = block_size*block_size;
    this->_gptr = matrix.begin().dart_gptr();
    this->_gptr.unitid = local_idx.unit.id;
    this->_gptr.addr_or_offs.offset = local_idx.index*sizeof(value_t);
#if 0
    if (!(local_idx.unit.id == (matrix.begin()+_glob_idx).dart_gptr().unitid)) {
      std::cout << "ERROR: block {" << block_row_idx << ","
                << block_col_idx << "} has two units: " << local_idx.unit.id
                << " vs " << (matrix.begin()+_glob_idx).dart_gptr().unitid << "!\n";
      std::cout << "gptr: " << (matrix.begin()+_glob_idx) << "\n";
      std::cout << "gptr: " << (matrix.begin()+_glob_idx).dart_gptr() << "\n";
      std::cout << "glob_idx: " << (_glob_idx) << "\n";
    }
#endif
    //assert(local_idx.unit.id == (matrix.begin()+_glob_idx).dart_gptr().unitid);
  }

  ~MatrixBlock() {
    if (!this->_is_local && _local_ptr != nullptr) {
#ifdef DEBUG
      std::cout << "Freeing pointer " << _local_ptr << std::endl;
#endif
      free(_local_ptr);
      _local_ptr = nullptr;
    }
  }

  MatrixBlock(const MatrixBlock& other) = default;

  MatrixBlock(MatrixBlock&& other) = delete;

  MatrixBlock&
  operator=(const MatrixBlock& other) = default;

  MatrixBlock&
  operator=(MatrixBlock&& other) = delete;


#if 0
  // TODO: this is currently not supported as `matrix->begin() + X` is broken...
  typename MatrixT::iterator
  begin() const {
    return _matrix->begin() + _glob_idx;
  }

  typename MatrixT::iterator
  end() const {
    return _matrix->begin() + _glob_idx + _size;
  }
#endif

  dart_gptr_t
  dart_gptr() const {
    return this->_gptr;
  }

  value_t *lbegin() {
    return this->local_ptr();
  }

  value_t *lend() {
    return this->local_ptr() + this->_size;
  }

  void fetch() {
    fetch_data();
  }

  bool is_local() {
    return this->_is_local;
  }

  size_t size() {
    return this->_size;
  }

  void
  fetch_async() {
    if (this->_local_ptr == nullptr) {
      this->_local_ptr = static_cast<value_t*>(malloc(this->_size * sizeof(value_t)));
#ifdef DEBUG_OUTPUT
      std::cout << "Fetching async " << _size << " elements into "
                << _local_ptr << std::endl;
#endif // DEBUG_OUTPUT
      dash::dart_storage<value_t> ds(_size);
      dart_get_handle(
        this->_local_ptr, this->_gptr, ds.nelem, ds.dtype, &_handle);
    }
  }

  void
  fetch_async(value_t *target) {
#ifdef DEBUG
    std::cout << "Fetching async " << _size << " elements into "
              << target << std::endl;
#endif
    dash::dart_storage<value_t> ds(_size);
    dart_get_handle(
      target, this->_gptr, ds.nelem, ds.dtype, ds.dtype, &_handle);
  }


  bool
  test()
  {
    if (_handle != DART_HANDLE_NULL) {
      int32_t flag;
      dart_test_local(&_handle, &flag);
      return flag != 0;
    }
    return true;
  }

  void wait()
  {
    if (_handle != DART_HANDLE_NULL) {
      dart_wait_local(&_handle);
    }
  }

  void store() {
    if (!_is_local && this->_local_ptr != NULL) {
      dash::dart_storage<value_t> ds(_size);
      dart_put_blocking(this->_gptr, _local_ptr, ds.nelem, ds.dtype, ds.dtype);
    }
  }

  void store_async() {
    if (!_is_local && this->_local_ptr != NULL) {
      dash::dart_storage<value_t> ds(_size);
      dart_put_handle(this->_gptr, _local_ptr, ds.nelem, ds.dtype, ds.dtype, &_handle);
    }
  }

  void store(value_t *source) {
    dash::dart_storage<value_t> ds(_size);
    dart_put_blocking(this->_gptr, source, ds.nelem, ds.dtype, ds.dtype);
  }

  void store_async(value_t *source) {
    if (!_is_local) {
      dash::dart_storage<value_t> ds(_size);
      dart_put_handle(this->_gptr, source, ds.nelem, ds.dtype, ds.dtype, &_handle);
    }
  }


  int unit() {
    return this->_matrix->pattern().unit_at(_glob_idx);
  }

  dart_handle_t dart_handle() {
    return this->_handle;
  }

private:

  void fetch_data() {
    if (!this->is_local()) {
      if (this->_local_ptr == nullptr) {
        this->_local_ptr = static_cast<value_t*>(malloc(this->_size * sizeof(value_t)));
      }
#ifdef DEBUG
      std::cout << "Fetching " << _size << " elements into "
                << _local_ptr << std::endl;
#endif
      dash::dart_storage<value_t> ds(_size);
      dart_get_blocking(
        this->_local_ptr, this->_gptr, ds.nelem, ds.dtype, ds.dtype);
        //dash::copy(begin, end, this->_local_ptr);
    }
  }

  value_t *local_ptr()
  {
    if (_is_local && _local_ptr == nullptr) {
      void *vptr;
      dart_gptr_getaddr(this->dart_gptr(), &vptr);
      _local_ptr = (value_t*)vptr;
      if (_local_ptr == nullptr) {
        std::cerr << "ERROR: failed to query local address for local block!" << std::endl;
      }
    }
    return _local_ptr;
  }

private:
  MatrixT* _matrix;
  value_t* _local_ptr = nullptr;
  int      _size = 0;
  int      _glob_idx;
  dart_handle_t _handle = DART_HANDLE_NULL;
  dart_gptr_t _gptr;
  bool     _is_local  = true;
};

#undef DEBUG

#endif // CHOLESKY_MATRIXBLOCK_H
