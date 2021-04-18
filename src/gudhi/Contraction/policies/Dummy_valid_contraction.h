/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_
#define CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_

#include <gudhi/Contraction/policies/Valid_contraction_policy.h>

namespace Gudhi {

namespace contraction {

/**
 *@brief Policy that accept all edge contraction.
 */
template< typename EdgeProfile>
class Dummy_valid_contraction : public Valid_contraction_policy<EdgeProfile> {
 public:
  typedef typename EdgeProfile::Point Point;

  bool operator()(const EdgeProfile& profile, const boost::optional<Point>& placement) {
    return true;
  }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_DUMMY_VALID_CONTRACTION_H_
