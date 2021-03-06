/*
 * Copyright 2018-2019 NWO-I
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

#include "Utility.h"

namespace Serialize {

  // Trivially copyable type: return sizeof
  template<typename T, typename std::enable_if<Detail::is_trivial<T>::value>::type* = nullptr>
  constexpr std::size_t serialize_size(T&&)
  {
    return sizeof(T);
  }

  // Ensure this is the worst match possible
  template<class T>
  constexpr std::size_t serialize_size(...)
  {
    return 0;
  }

  // Version for string
  template<class C>
  std::size_t serialize_size(const std::basic_string<C>& t)
  {
    return t.size() * sizeof(C);
  }

} // namespace Serialize

namespace Detail {

  // Specialization for std::tuple
  template<std::size_t I, class... Args>
  struct tuple_size_helper {
    static constexpr std::size_t impl(const std::tuple<Args...>& t)
    {
      return Serialize::serialize_size(std::get<I - 1>(t)) + tuple_size_helper<I - 1, Args...>::impl(t);
    }
  };

  template<class... Args>
  struct tuple_size_helper<0, Args...> {
    static constexpr std::size_t impl(const std::tuple<Args...>&) { return 0; }
  };

} // namespace Detail

namespace Serialize {

  // Version for tuple
  template<class... Args>
  constexpr std::size_t serialize_size(const std::tuple<Args...>& t)
  {
    return Detail::tuple_size_helper<sizeof...(Args), Args...>::impl(t);
  }

  // Version for pair
  template<class P, typename std::enable_if<is_pair<P>::value>::type* = nullptr>
  std::size_t serialize_size(const P& p)
  {
    return serialize_size(p.first) + serialize_size(p.second);
  }

  // Version for container
  template<
    class C,
    typename std::enable_if<is_container<C>::value && Detail::is_trivial<typename C::value_type>::value>::type* =
      nullptr>
  std::size_t serialize_size(const C& c)
  {
    return sizeof(typename C::value_type) * c.size();
  }

  template<
    class C,
    typename std::enable_if<is_container<C>::value && !Detail::is_trivial<typename C::value_type>::value>::type* =
      nullptr>
  std::size_t serialize_size(const C& c)
  {
    std::size_t s = 0;
    for (const auto& e : c) {
      s += serialize_size(e);
    }
    return s;
  }

  // Specialization for std::array
  template<class T, std::size_t N, typename std::enable_if<Detail::is_trivial<T>::value>::type* = nullptr>
  constexpr std::size_t serialize_size(const std::array<T, N>& v)
  {
    return sizeof(T) * N;
  }

  template<class T, std::size_t N, typename std::enable_if<!Detail::is_trivial<T>::value>::type* = nullptr>
  std::size_t serialize_size(const std::array<T, N>& v)
  {
    std::size_t s = 0;
    for (const auto& e : v) {
      s += serialize_size(e);
    }
    return s;
  }

} // namespace Serialize
