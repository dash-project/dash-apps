#ifndef UNITIALIZED_VECTOR_H_
#define UNITIALIZED_VECTOR_H_

#include <utility>
#include <cstdlib>

template<typename T>
struct uninitialized_vector {
public:

  typedef uninitialized_vector<T> self_t;

  uninitialized_vector(size_t size) : _size(size) {
    this->_data = static_cast<T*>(malloc(size*sizeof(T)));
  }

  uninitialized_vector(self_t&& other) noexcept
    : _data(other._data), _size(other._size) {
    other._data = NULL;
    other._size = 0;
  }

  // no, we don't allow copying (for now)
  uninitialized_vector(const self_t& other) = delete;

  ~uninitialized_vector() {
    this->free();
  }

  self_t&
  operator=(self_t&& other) noexcept {
    std::swap(this->_data, other._data);
    std::swap(this->_size, other._size);
    return *this;
  }

  // no, we don't allow copying (for now)
  self_t&
  operator=(const self_t& other) = delete;

  T& operator[](size_t pos) noexcept {
    return this->_data[pos];
  }

  const T& operator[](size_t pos) const noexcept {
    return this->_data[pos];
  }


  T* data(void) noexcept {
    return this->_data;
  }

  const T* data(void) const noexcept {
    return this->_data;
  }

  void free() {
    ::free(this->_data);
    this->_data = NULL;
    this->_size = 0;
  }

        T* begin()       noexcept { return data(); }
        T* end()         noexcept { return data() + _size; }

  const T* begin() const noexcept { return data(); }
  const T* end()   const noexcept { return data() + _size; }

private:
  T *_data;
  size_t _size;
};



#endif /* UNITIALIZED_VECTOR_H_ */
