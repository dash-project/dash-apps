#ifndef CHOLESKY_MATRIXBLOCK_H
#define CHOLESKY_MATRIXBLOCK_H

#include <libdash.h>

template<typename MatrixT>
class MatrixBlock {
public:

  using value_t = typename MatrixT::value_type;

  MatrixBlock(
    MatrixT &matrix,
    size_t block_row_idx,
    size_t block_col_idx)
  : _matrix(&matrix) {
    auto& pattern       = matrix.pattern();
    auto  glob_index    = block_row_idx*pattern.blocksize(0)*pattern.extent(1) + block_col_idx*pattern.blocksize(1);
    this->_glob_idx     = glob_index;
    this->_is_local = pattern.is_local(glob_index);
    if (_is_local) {
      _local_ptr = matrix.lbegin() +
                    pattern.local_index(
                      {block_row_idx*pattern.blocksize(0),
                        block_col_idx*pattern.blocksize(1)}).index;
    }
    this->_size = pattern.blocksize(0) * pattern.blocksize(1);
#ifdef DEBUG
    std::cout << "BlockCache: " << block_row_idx << "x" << block_col_idx << " "
              << pattern.blocksize(0) << "x" << pattern.blocksize(1)
              << " (size=" << _size << ", is_local=" << _is_local
              << ", glob_idx=" << glob_index
              << std::endl;
#endif
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

  MatrixBlock(const MatrixBlock& other)
  : _matrix(other._matrix),
    _local_ptr(other._is_local ? other._local_ptr : nullptr),
    _size(other._size),
    _glob_idx(other._glob_idx),
    _is_local(other._is_local)
  { }

  MatrixBlock(MatrixBlock&& other) = delete;

  MatrixBlock&
  operator=(const MatrixBlock& other)
  {
    if (this == &other) return *this;

    _matrix = other._matrix;
    _local_ptr = other._is_local ? other._local_ptr : nullptr;
    _size      = other._size;
    _glob_idx  = other._glob_idx;
    _is_local  = other._is_local;
  }

  MatrixBlock&
  operator=(MatrixBlock&& other) = delete;


  typename MatrixT::iterator
  begin() const {
    return _matrix->begin() + _glob_idx;
  }

  typename MatrixT::iterator
  end() const {
    return _matrix->begin() + _glob_idx + _size;
  }

  value_t *lbegin() {
    if (!this->_is_local) {
      fetch_data();
    }
    return this->_local_ptr;
  }

  value_t *lend() {
    if (!this->_is_local) {
      fetch_data();
      return this->_local_ptr + this->_size;
    }
    return this->_local_ptr + this->_size;
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
      auto begin = this->_matrix->begin() + _glob_idx;
#ifdef DEBUG
      std::cout << "Fetching async " << _size << " elements into "
                << _local_ptr << std::endl;
#endif
      dash::dart_storage<value_t> ds(_size);
      dart_get_handle(
        this->_local_ptr, begin.dart_gptr(), ds.nelem, ds.dtype, &_handle);
      //dash::copy(begin, end, this->_local_ptr);
    }
  }

  void
  fetch_async(value_t *target) {
    auto begin = this->_matrix->begin() + _glob_idx;
#ifdef DEBUG
    std::cout << "Fetching async " << _size << " elements into "
              << target << std::endl;
#endif
    dash::dart_storage<value_t> ds(_size);
    dart_get_handle(
      target, begin.dart_gptr(), ds.nelem, ds.dtype, &_handle);
    //dash::copy(begin, end, this->_local_ptr);
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
      auto begin = this->_matrix->begin() + _glob_idx;
      dash::dart_storage<value_t> ds(_size);
      dart_put_blocking(begin.dart_gptr(), _local_ptr, ds.nelem, ds.dtype);
    }
  }

private:

  void fetch_data() {
    if (this->_local_ptr == nullptr) {
      this->_local_ptr = static_cast<value_t*>(malloc(this->_size * sizeof(value_t)));
      auto begin = this->_matrix->begin() + _glob_idx;
#ifdef DEBUG
      std::cout << "Fetching " << _size << " elements into "
                << _local_ptr << std::endl;
#endif
      dash::dart_storage<value_t> ds(_size);
      dart_get_blocking(
        this->_local_ptr, begin.dart_gptr(), ds.nelem, ds.dtype);
      //dash::copy(begin, end, this->_local_ptr);
    }
  }

private:
  MatrixT* _matrix;
  value_t* _local_ptr = nullptr;
  size_t   _size = 0;
  size_t   _glob_idx;
  dart_handle_t _handle = DART_HANDLE_NULL;
  bool     _is_local  = true;
};

#endif // CHOLESKY_MATRIXBLOCK_H
