/*
 * src/tracker/analysis/vertex.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <tracker/analysis/digi.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "../helper/analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////



//__Vertex Default Constructor__________________________________________________________________
digi::digi() {
  clear();
}
//----------------------------------------------------------------------------------------------

//__Vertex Constructor__________________________________________________________________________
digi::digi(const analysis::digi_event& points) {
  reset(points);
}
//----------------------------------------------------------------------------------------------




//__Clear Hits from Track_______________________________________________________________________
void digi::clear() {
  reset(analysis::digi_event{});
}
//----------------------------------------------------------------------------------------------




//__Fill Plots with Vertexing Variables_________________________________________________________
void digi::fill_plots(plot::histogram_collection& collection,
                        const digi::plotting_keys& keys) const {
  if (collection.count(keys.t)) collection[keys.t].insert(digi_t_value() / units::time);
  if (collection.count(keys.x)) collection[keys.x].insert(digi_x_value() / units::length);
  if (collection.count(keys.y)) collection[keys.y].insert(digi_y_value() / units::length);
  if (collection.count(keys.z)) collection[keys.z].insert(digi_z_value() / units::length);
  if (collection.count(keys.e)) collection[keys.e].insert(digi_e_value() / units::energy);

  if (collection.count(keys.px)) collection[keys.px].insert(digi_px_value() / units::momentum);
  if (collection.count(keys.py)) collection[keys.py].insert(digi_py_value() / units::momentum);
  if (collection.count(keys.pz)) collection[keys.pz].insert(digi_pz_value() / units::momentum);

}
//----------------------------------------------------------------------------------------------




//__Vertex Data Tree Constructor________________________________________________________________
digi::tree::tree(const std::string& name)
    : tree(name, name) {}
//----------------------------------------------------------------------------------------------

//__Vertex Data Tree Constructor________________________________________________________________
digi::tree::tree(const std::string& name,
                   const std::string& title)
    : analysis::tree(name, title),
      t(emplace_branch<real_branch_value_type>("Digi_t")),
      x(emplace_branch<real_branch_value_type>("Digi_x")),
      y(emplace_branch<real_branch_value_type>("Digi_y")),
      z(emplace_branch<real_branch_value_type>("Digi_z")),
      e(emplace_branch<real_branch_value_type>("Digi_energy")),
      px(emplace_branch<real_branch_value_type>("Digi_px")),
      py(emplace_branch<real_branch_value_type>("Digi_py")),
      pz(emplace_branch<real_branch_value_type>("Digi_pz")),
      _vector_branches({t, x, y, z, e, px, py, pz}) {}
//----------------------------------------------------------------------------------------------

//__Track Data Tree Insertion___________________________________________________________________
void digi::tree::insert(const digi& digi) {
  t.get().push_back(digi.digi_t_value() / units::time);
  x.get().push_back(digi.digi_x_value() / units::length);
  y.get().push_back(digi.digi_y_value() / units::length);
  z.get().push_back(digi.digi_z_value() / units::length);
  e.get().push_back(digi.digi_e_value() / units::energy);
  px.get().push_back(digi.digi_px_value() / units::momentum);
  py.get().push_back(digi.digi_py_value() / units::momentum);
  pz.get().push_back(digi.digi_pz_value() / units::momentum);

}
//----------------------------------------------------------------------------------------------

//__Clear Vertex Data Tree______________________________________________________________________
void digi::tree::clear() {
  _count = 0ULL;
  for (auto& entry : _vector_branches)
    entry.get().get().clear();

}
//----------------------------------------------------------------------------------------------

//__Reserve Space for Vertex Data Tree__________________________________________________________
void digi::tree::reserve(std::size_t capacity) {
  for (auto& entry : _vector_branches)
    entry.get().get().reserve(capacity);
}
//----------------------------------------------------------------------------------------------



} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
