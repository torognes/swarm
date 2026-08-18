/*
  SWARM

  Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#ifndef SWARM_UTILS_MAKE_UNIQUE_H
#define SWARM_UTILS_MAKE_UNIQUE_H


#include <memory>   // std::unique_ptr
#include <utility>  // std::forward


namespace utils {

  // returns a std::unique_ptr for a given type (following http://herbsutter.com/gotw/_102/)
  // always call as: utils::make_unique<>() to avoid name overloads if compiled with C++17.
 
  template<typename Type, typename... Args>
  auto make_unique(Args&&... args) -> std::unique_ptr<Type> {
    return std::unique_ptr<Type>(new Type(std::forward<Args>(args)...));
  }

} // namespace utils

#endif  // SWARM_UTILS_MAKE_UNIQUE_H
