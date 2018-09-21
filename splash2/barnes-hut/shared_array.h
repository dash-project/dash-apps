#ifndef SHARED_ARRAY_H__INCLUDED
#define SHARED_ARRAY_H__INCLUDED

#include <dash/Array.h>
#include <dash/algorithm/Copy.h>

namespace dash {

template <class T>
class SharedArray {
  static_assert(std::is_array<T>::value, "only allowed for array types");
  static_assert(
      std::rank<T>::value < 2, "only 1 dimensional arrays are supported");

  using underlying_value_t = typename std::remove_all_extents<T>::type;

  using pattern_t = dash::CSRPattern<1>;
  using container_t   = dash::
      Array<underlying_value_t, typename pattern_t::index_type, pattern_t>;
  using size_type = pattern_t::size_type;

  static constexpr const size_t extent = std::extent<T>::value;

public:
  using value_type      = T;
  using reference       = dash::GlobRef<underlying_value_t>;
  using const_reference = typename reference::const_type;

public:
  SharedArray(
      team_unit_t owner = team_unit_t{0},
      Team const& team  = dash::Team::All())
    : m_team(&team)
    , m_owner(owner)
    , m_container(const_cast<dash::Team&>(team))
  {
    if (dash::is_initialized()) {
      DASH_ASSERT_RETURNS(init(), true);
    }
  }

  bool init()
  {
    std::vector<typename pattern_t::size_type> sizes;
    sizes.reserve(m_team->size());

    std::fill_n(std::back_inserter(sizes), m_owner, 0);
    sizes.push_back(extent);
    std::fill_n(std::back_inserter(sizes), m_team->size() - m_owner - 1, 0);

    pattern_t pattern{sizes};

    DASH_ASSERT_EQ(pattern.size(), m_team->size(), "invalid pattern size");

    return m_container.allocate(pattern);
  }

  reference operator[](size_type global_index)
  {
    DASH_ASSERT(global_index < extent);
    return m_container[global_index];
  }

  reference at(size_type global_index)
  {
    if (global_index >= extent) {
      throw std::out_of_range("out of bounds access in shared_array");
    }
    return m_container[global_index];
  }

  const_reference operator[](size_type global_index) const
  {
    DASH_ASSERT(global_index < extent);
    return m_container[global_index];
  }

  const_reference at(size_type global_index) const
  {
    if (global_index >= extent) {
      throw std::out_of_range("out of bounds access in shared_array");
    }
    return m_container[global_index];
  }

  void set(value_type const& val)
  {
    if (m_container.begin().is_local()) {
      // local copy
      std::copy(std::begin(val), std::end(val), m_container.lbegin());
    }
    else {
      // global copy
      std::vector<dart_handle_t> handles;
      dash::internal::copy_impl(
          &(*std::begin(val)), &(*std::end(val)), m_container.begin(), handles);

      if (!handles.empty()) {
        dart_waitall(handles.data(), handles.size());
      }
    }
  }

  void get(value_type& out) const
  {
    if (m_container.begin().is_local()) {
      std::copy(m_container.lbegin(), m_container.lend(), std::begin(out));
    }
    else {
      std::vector<dart_handle_t> handles;

      dash::internal::copy_impl(
          m_container.begin(),
          m_container.begin() + extent,
          &(*std::begin(out)),
          handles);

      if (!handles.empty()) {
        dart_waitall(handles.data(), handles.size());
      }
    }
  }

private:
  dash::Team const* m_team;
  dash::team_unit_t m_owner{0};
  container_t           m_container;
};


//Definition of static member
template<typename T>
constexpr size_t SharedArray<T>::extent;

}  // namespace dash
#endif
