/*
 * include/tracker/analysis/vertex.hh
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

#ifndef TRACKER__ANALYSIS__DIGI_HH
#define TRACKER__ANALYSIS__DIGI_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/analysis/tree.hh>
#include <tracker/geometry.hh>
#include <tracker/plot.hh>


namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Vertex Object_______________________________________________________________________________
class digi {
public:
  class tree;
  enum class parameter { T, X, Y, Z, E, PX, PY, PZ };


  using container_type = analysis::digi_event;
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;
  using reverse_iterator = typename container_type::reverse_iterator;
  using const_reverse_iterator = typename container_type::const_reverse_iterator;


  digi();
  digi(const digi_event& points);
  digi(digi_event&& points);
  digi(const digi& rhs)     = default;
  digi(digi&& rhs) noexcept = default;
  digi& operator=(const digi& rhs)     = default;
  digi& operator=(digi&& rhs) noexcept = default;





  const digi_hit point() const;



  real digi_t_value()   const { return _digi_hit.t; }
  real digi_x_value()  const { return _digi_hit.x; }
  real digi_y_value()  const { return _digi_hit.y; }
  real digi_z_value()  const { return _digi_hit.z; }
  real digi_e_value()  const { return _digi_hit.e; }
  real digi_px_value() const { return _digi_hit.px; }
  real digi_py_value() const { return _digi_hit.py; }
  real digi_pz_value() const { return _digi_hit.pz; }
  real value(const parameter p) const;





  // const analysis::digi_event& digi_event() const noexcept { return _digi_event; }

  std::size_t size() const { return _digi_event.size(); }
  bool empty() const noexcept { return _digi_event.empty(); }


  iterator       begin()        noexcept { return _digi_event.begin();  }
  const_iterator begin()  const noexcept { return _digi_event.cbegin(); }
  iterator       end()          noexcept { return _digi_event.end();    }
  const_iterator end()    const noexcept { return _digi_event.cend();   }
  const_iterator cbegin() const noexcept { return _digi_event.cbegin(); }
  const_iterator cend()   const noexcept { return _digi_event.cend();   }

  reverse_iterator       rbegin()        noexcept { return _digi_event.rbegin();  }
  const_reverse_iterator rbegin()  const noexcept { return _digi_event.crbegin(); }
  reverse_iterator       rend()          noexcept { return _digi_event.rend();    }
  const_reverse_iterator rend()    const noexcept { return _digi_event.crend();   }
  const_reverse_iterator crbegin() const noexcept { return _digi_event.crbegin(); }
  const_reverse_iterator crend()   const noexcept { return _digi_event.crend();   }



  std::size_t reset(const analysis::digi_event& points);

  std::size_t insert(const digi_hit& point);
  std::size_t insert(const analysis::digi_event& points);

  std::size_t remove(const std::size_t index);
  std::size_t remove(const std::vector<std::size_t>& indices);

  void clear();

  struct plotting_keys {
    plot::histogram::name_type t, x, y, z, e, px, py, pz;
  };

  void fill_plots(plot::histogram_collection& collection,
                  const plotting_keys& keys) const;

  void draw(plot::canvas& canvas,
            const real size,
            const plot::color color,
            const bool with_errors=false) const;


  bool operator==(const digi& other) const noexcept { return _digi_event == other._digi_event; }
  bool operator!=(const digi& other) const noexcept { return !(*this == other);        }



protected:
  analysis::digi_hit _digi_hit;
  analysis::digi_event _digi_event;
};
//----------------------------------------------------------------------------------------------

//__Vector of digi Hits__________________________________________________________________________
using digi_vector = std::vector<digi>;
//----------------------------------------------------------------------------------------------

//__Track Data Tree Specialization______________________________________________________________
class digi::tree : public analysis::tree {
public:
  using real_branch_value_type = std::vector<double>;
  using real_branch_type = branch<real_branch_value_type>;

  tree(const std::string& name);
  tree(const std::string& name,
       const std::string& title);

  real_branch_type t, x, y, z, e, px, py, pz;


  void insert(const digi& digi);
  void clear();
  void reserve(std::size_t capacity);


  template<class UnaryPredicate>
  UnaryPredicate fill_if(const digi_vector& digis,
                         UnaryPredicate f) {
    clear();
    reserve(digis.size());
    for (const auto& digi : digis) {
      if (f(digi))
        insert(digi);
    }
    analysis::tree::fill();
    return std::move(f);
  }

  void fill(const digi_vector& digis={}) {
    fill_if(digis, [](auto) { return true; });
  }



private:
  branch<unsigned long long> _count;
  std::vector<std::reference_wrapper<real_branch_type>> _vector_branches;
};
//----------------------------------------------------------------------------------------------



} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__VERTEX_HH */
