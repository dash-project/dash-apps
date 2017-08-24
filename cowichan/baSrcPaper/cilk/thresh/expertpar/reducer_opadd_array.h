#ifndef __REDUCER_OPADD_ARRAY_H_
#define __REDUCER_OPADD_ARRAY_H_

#include <cilk/reducer_opadd.h>

// An object which encapsulates an array T[SIZE];
template <typename T, int SIZE>
  class reducer_opadd_array_view {
 public:
  inline T& operator[](int idx) {
    return m_a[idx];
  }

  // Try for 128-byte alignment and padding.
  __declspec(align(128)) T m_a[SIZE];
  char PAD[128];
};

template <typename T, int SIZE>
  class reducer_opadd_array {
 public:
  class Monoid : public cilk::monoid_base<reducer_opadd_array_view<T, SIZE> > {
  public:
    static void reduce(reducer_opadd_array_view<T, SIZE>* left,
                       reducer_opadd_array_view<T, SIZE>* right) {
      left->m_a[0:SIZE] += right->m_a[0:SIZE];
    }
  };

  // Get the view of the object (for parallel access).
  T* get_view(void) {
    return imp_.view().m_a;
  }

  T& operator[](int idx) {
    return imp_.view()[idx];
  }
    
 private:
  // Hyperobject for views.
  cilk::reducer<Monoid> imp_;
};
#endif
