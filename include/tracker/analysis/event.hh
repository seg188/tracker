/*
 * include/tracker/analysis/event.hh
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

#ifndef TRACKER__ANALYSIS__EVENT_HH
#define TRACKER__ANALYSIS__EVENT_HH
#pragma once

#include <tracker/analysis/type.hh>
#include <tracker/geometry.hh>
#include <tracker/util/algorithm.hh>
#include <iostream>
#include <vector>

//__For Printing Type std::vector_______________________________________________________________
namespace std {
	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const vector<T>& v)
	{
		os << "[";
		for (int i = 0; i < v.size(); ++i) {
			os << v[i];
			if (i != v.size() - 1)
				os << ", ";
		}
		os << "]\n";
		return os;
	}
}
//----------------------------------------------------------------------------------------------

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Reduce Event Vector to Event________________________________________________________________
const event reduce(const event_vector& events);
const full_event reduce(const full_event_vector& full_events);
const energy_event reduce(const energy_event_vector& energy_events);
const complete_event reduce(const complete_event_vector& complete_events);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per unit Length____________________________________________________
const r4_point event_density(const event& points);
const r4_point event_density(const full_event& points);
const r4_point event_density(const energy_event& energy_points);
const r4_point event_density(const complete_event& complete_points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per unit Volume____________________________________________________
// TODO: real event_volume_density(const event& points);
// TODO: real event_volume_density(const full_event& points);
//----------------------------------------------------------------------------------------------

//__Calculate Number of Hits per Geometric Element______________________________________________
// TODO: real geometric_event_density(const event& points);
// TODO: real geometric_event_density(const full_event& points);
//----------------------------------------------------------------------------------------------

template<class Geometry=void>
const event add_digi_event(const event& points, const energy_event& energy_points, const complete_event& complete_points) {

  complete_event c_out = complete_points;
  complete_event c_time_sorted = t_sort(c_out);
  //std::cout << "CCCCCCCCCCCCCCCCCCCC: " << c_time_sorted;

  std::vector<long double> detector_ids;
  detector_ids.reserve(c_time_sorted.size());
  for (const auto& h : c_time_sorted) {
      detector_ids.push_back(h.det_id);
  }
  util::algorithm::sort_range(detector_ids);
  detector_ids.erase(std::unique(detector_ids.begin(), detector_ids.end()), detector_ids.cend());

  event digi_out;
  digi_out.reserve(detector_ids.size());


  const auto spacing = 20 * units::time;
  const auto threshold = 0.65 * units::energy;

  for (const auto& d : detector_ids) {
	  long double energy_sum = 0;
      long double weighted_time = 0;
      long double weighted_z = 0;
	  long double point_x = 0;
      long double point_y = 0;

	  int counter = 0;
	  std::vector<long double> times;
	  times.reserve(counter);
	  std::vector<long double> zs;
      zs.reserve(counter);
	  std::vector<long double> ys;
      ys.reserve(counter);
	  std::vector<long double> xs;
      xs.reserve(counter);
	  std::vector<long double> deposits;
      deposits.reserve(counter);

      for (const auto& h : c_time_sorted) {
          if (h.det_id == d){
              times.push_back(h.t);
              zs.push_back(h.z);
              xs.push_back(h.x);
              ys.push_back(h.y);
              deposits.push_back(h.e);
			  ++counter;
		  }
	  }

	  for (const auto& h : c_time_sorted) {
		  if (h.det_id == d){
			  auto starting_index = 0;
			  while (starting_index < times.size()) {
				  auto t0 = times[starting_index];
				  std::vector<long double> e_components;
				  std::vector<long double> t_components;
				  std::vector<long double> z_components;
				  for (int i = 0; i < times.size(); i++) {
					  if (times[i] < (t0 + spacing)) {
						  e_components.push_back(deposits[i]);
						  t_components.push_back(times[i]*deposits[i]);
                          z_components.push_back(zs[i]*deposits[i]);
					  }
				  }
				  long double e_sum = std::accumulate(e_components.begin(), e_components.end(), 0.0L);
				  long double t_sum = std::accumulate(t_components.begin(), t_components.end(), 0.0L);
				  long double z_sum = std::accumulate(z_components.begin(), z_components.end(), 0.0L);
				  if (e_sum > threshold) {
                      energy_sum = e_sum;
					  weighted_time = t_sum;
                      weighted_z = z_sum;
                      starting_index += e_components.size();
				  } else {
                      starting_index += 1;
				  }
			  }

			  point_x = xs[0];
			  point_y = ys[0];

		  }
	  }
	  //need to divide by 10 for some reason!
	  digi_out.push_back({weighted_time/energy_sum * units::time, (point_x/10) * units::length, (point_y/10) * units::length, ((weighted_z/energy_sum)/10) * units::length});

  }
  //  std::cout << "OOOOOOOOOOOOOOOOOOOOOOO: " << out <<std::endl;
  //  std::cout << "NNNNNNNNNNNNNNNNNNNNNNN: " << digi_out <<std::endl;
  return digi_out;
}


// template<class Geometry=void>
// const full_hit add_digi(const hit& point, const energy_hit& energy_point) {
//   const auto volume = geometry::custom::volume<Geometry>(point);
//   const auto limits = geometry::custom::limits_of<Geometry>(volume);
//   const auto center = limits.center;
//   const auto min = limits.min;
//   const auto max = limits.max;
//   return full_hit{
// 	point.t, center.x, center.y, point.z,
// 		{geometry::custom::time_resolution_of<Geometry>(volume),
// 				max.x - min.x,
// 				max.y - min.y,
// 				max.z - min.z}};
// }
// template<class Geometry=void>
// const full_event add_digi(const event& points, const energy_event& energy_points) {
//   full_event out;
//   out.reserve(points.size());
//   util::algorithm::back_insert_transform_two(points, energy_points, out,
//     [](const auto& point, const auto& energy_point) { return add_digi<Geometry>(point, energy_point); });
//   return out;
// }

//__Find The Errors Associated with a Hit from Geometry_________________________________________
template<class Geometry=void>
const full_hit add_width(const hit& point) {
  const auto volume = geometry::custom::volume<Geometry>(point);
  const auto limits = geometry::custom::limits_of<Geometry>(volume);
  const auto center = limits.center;
  const auto min = limits.min;
  const auto max = limits.max;
  return full_hit{
    point.t, center.x, center.y, point.z,
    {geometry::custom::time_resolution_of<Geometry>(volume),
     max.x - min.x,
     max.y - min.y,
     max.z - min.z}};
}
template<class Geometry=void>
const full_event add_width(const event& points) {
  full_event out;
  out.reserve(points.size());
  util::algorithm::back_insert_transform(points, out,
    [](const auto& point) { return add_width<Geometry>(point); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Event Partition Type________________________________________________________________________
struct event_partition {
  event_vector parts; Coordinate coordinate; real interval;
};
struct full_event_partition {
  full_event_vector parts; Coordinate coordinate; real interval;
};
struct energy_event_partition {
	energy_event_vector parts; Coordinate coordinate; real interval;
};
//----------------------------------------------------------------------------------------------

//__Partition Points by Coordinate______________________________________________________________
const event_partition partition(const event& points,
                                const Coordinate coordinate,
                                const real interval);
const full_event_partition partition(const full_event& points,
                                     const Coordinate coordinate,
                                     const real interval);
const energy_event_partition partition(const energy_event& points,
									 const Coordinate coordinate,
									 const real interval);
//----------------------------------------------------------------------------------------------

//__Reset Partition by new Interval_____________________________________________________________
const event_partition repartition(const event_partition& previous,
                                  const real interval);
const full_event_partition repartition(const full_event_partition& previous,
                                       const real interval);
const energy_event_partition repartition(const energy_event_partition& previous,
										   const real interval);
//----------------------------------------------------------------------------------------------

//__Reduce Event Partition to Events____________________________________________________________
const event reduce_partition(const event_partition& previous);
const full_event reduce_partition(const full_event_partition& previous);
const energy_event reduce_partition(const energy_event_partition& previous);
//----------------------------------------------------------------------------------------------

//__Calculate Density of Partition______________________________________________________________
// TODO: real partition_density(const event_partition& points);
// TODO: real partition_density(const full_event_partition& points);
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__EVENT_HH */
