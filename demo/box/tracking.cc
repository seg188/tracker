/*
 * demo/box/tracking.cc
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

#include <tracker/analysis.hh>

#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>
#include <tracker/script.hh>

#include <tracker/analysis/type.hh>
#include <tracker/core/type.hh>
#include <tracker/core/units.hh>

#include "geometry.hh"
#include "io.hh"
#include "TTree.h"
#include <TFile.h>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class std::vector<double>+;
#endif


//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace mc       = analysis::mc;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
namespace script   = MATHUSLA::TRACKER::script;
namespace reader   = MATHUSLA::TRACKER::reader;
namespace io       = MATHUSLA::box::io;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {



//__Alter Event_________________________________________________________________________________
const analysis::full_event alter_event(const analysis::full_event& event,
                                       const script::tracking_options& options,
                                       const box::io::extension_parser& extension) {
  return box::geometry::restrict_layer_count(
           mc::use_efficiency(
             mc::add_noise<box::geometry>(
               event,
               options.simulated_noise_rate,
               options.event_time_window.begin,
               options.event_time_window.end,
               box::geometry::full(extension.layer_count)),
             options.simulated_efficiency),
           extension.layer_count);
}
//----------------------------------------------------------------------------------------------

//__Find Tracks for Box_________________________________________________________________________
const analysis::track_vector find_tracks(const analysis::full_event& event,
                                         const script::tracking_options& options,
                                         const type::real limit_chi_squared,
                                         const std::size_t overlap,
                                         plot::canvas& canvas,
                                         analysis::full_event& non_track_points) {
  namespace ash = analysis::seed_heuristic;
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(
    options.seed_size,
    layers,
    ash::all<ash::double_cone, ash::piecewise_speed>{
      {options.line_width},
      {0.9L * units::speed_of_light}});

  /*
  for (const auto& seed : seeds) {
    for (std::size_t i{}; i < seed.size() - 1UL; ++i)
      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i + 1UL]), 1, plot::color::BLACK);
    for (const auto& point : seed)
      std::cout << type::reduce_to_r4(point) << " ";
    std::cout << "\n";
  }
  */

  const auto joined = analysis::join_all(seeds);
  auto first_tracks = analysis::independent_fit_seeds(joined, options.layer_axis);

  for (auto& track : first_tracks)
    track.prune_on_chi_squared(limit_chi_squared);

  const auto out = analysis::overlap_fit_tracks(first_tracks, overlap);
  non_track_points = analysis::non_tracked_points(event, out, true);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Track an Event Bundle_______________________________________________________________________
void track_event_bundle(const script::path_vector& paths,
                        const mc::full_event_vector_bundle& bundle,
                        const script::tracking_options& options,
                        const box::io::extension_parser& extension,
                        const script::path_type& save_path) {
  static const plot::value_tag filetype_tag("FILETYPE", "MATHUSLA TRACKING STATFILE");
  static const plot::value_tag project_tag("PROJECT", "Box");

  const auto energy_events = bundle.energy_events;
  const auto complete_events = bundle.complete_events;
  const auto imported_events = bundle.events;

  const auto import_size = imported_events.size();
  if (import_size == 0UL)
    return;

  const auto mc_imported_events = bundle.true_events;

  std::cout << "Event Count: " << import_size << "\n";

  analysis::track::tree track_tree{"track_tree", "MATHUSLA Track Tree"};
  analysis::vertex::tree vertex_tree{"vertex_tree", "MATHUSLA Vertex Tree"};

  // track_tree.set_file(save_path);
  // vertex_tree.set_file(save_path);
  track_tree.add_friend(vertex_tree, "vertex");
  vertex_tree.add_friend(track_tree, "track");

  //___create digi_tree_________________________________________________________________________
  TTree digi_tree("digi_tree", "MATHUSLA Digi Tree");

    std::vector<double> digi_hit_t;
    std::vector<double> digi_hit_x;
    std::vector<double> digi_hit_y;
    std::vector<double> digi_hit_z;
    std::vector<double> digi_hit_e;
    std::vector<double> digi_hit_px;
    std::vector<double> digi_hit_py;
    std::vector<double> digi_hit_pz;
    double Digi_numHits;

    auto branch_digi_num_hits  = digi_tree.Branch("Digi_numHits", &Digi_numHits, "Digi_numHits/D");
    auto branch_t  = digi_tree.Branch("Digi_time", "std::vector<double>", &digi_hit_t, 32000, 99);
    auto branch_x  = digi_tree.Branch("Digi_x", "std::vector<double>", &digi_hit_x, 32000, 99);
    auto branch_y  = digi_tree.Branch("Digi_y", "std::vector<double>", &digi_hit_y, 32000, 99);
    auto branch_z  = digi_tree.Branch("Digi_z", "std::vector<double>", &digi_hit_z, 32000, 99);
    auto branch_e  = digi_tree.Branch("Digi_energy", "std::vector<double>", &digi_hit_e, 32000, 99);
    auto branch_px = digi_tree.Branch("Digi_px", "std::vector<double>", &digi_hit_px, 32000, 99);
    auto branch_py = digi_tree.Branch("Digi_py", "std::vector<double>", &digi_hit_py, 32000, 99);
    auto branch_pz = digi_tree.Branch("Digi_pz", "std::vector<double>", &digi_hit_pz, 32000, 99);
  //____________________________________________________________________________________________



  for (std::size_t event_counter{}; event_counter < import_size; ++event_counter) {


    const auto digitized_full_event = analysis::full_digi_event<box::geometry>(imported_events[event_counter], energy_events[event_counter], complete_events[event_counter]);
    for (const auto& h : digitized_full_event) {
		digi_hit_e.push_back(h.e);
		digi_hit_px.push_back(h.px);
		digi_hit_py.push_back(h.py);
		digi_hit_pz.push_back(h.pz);
	}


    const auto digi_event = analysis::add_digi_event<box::geometry>(imported_events[event_counter], energy_events[event_counter], complete_events[event_counter]);

	const auto event = analysis::add_width<box::geometry>(digi_event);
	const auto event_size = event.size();
	const auto event_counter_string = std::to_string(event_counter);

	const auto compressed_event_t = options.time_smearing ? mc::time_smear<box::geometry>(mc::compress<box::geometry>(event))
	   	                                                  : mc::compress<box::geometry>(event);

	const auto compressed_event = options.positionz_smearing ? mc::positionz_smear<box::geometry>(compressed_event_t)
		                                                     : compressed_event_t;

    for (const auto& h : compressed_event) {
        digi_hit_t.push_back(h.t);
        digi_hit_x.push_back(h.x);
        digi_hit_y.push_back(h.y);
        digi_hit_z.push_back(h.z);
    }

    Digi_numHits = digi_hit_x.size();
    if (digitized_full_event.size() == 0UL && digi_event.size() == 0UL) {
      digi_tree.Fill();
      continue;
    }
    if (digitized_full_event.size() != 0UL && digi_event.size() != 0UL) {
      digi_tree.Fill();
    }
    digi_hit_t.clear();
    digi_hit_x.clear();
    digi_hit_y.clear();
    digi_hit_z.clear();
    digi_hit_e.clear();
    digi_hit_px.clear();
    digi_hit_py.clear();
    digi_hit_pz.clear();


    const auto compression_size = event_size / static_cast<type::real>(compressed_event.size());

    if (event_size == 0UL || compression_size == event_size) {
      track_tree.fill();
      vertex_tree.fill();
      continue;
    }

    const auto altered_event = alter_event(compressed_event, options, extension);
    const auto event_density = box::geometry::event_density(altered_event);
    box::io::print_event_summary(event_counter, altered_event.size(), compression_size, event_density);

    if (event_density >= options.event_density_limit) {
      track_tree.fill();
      vertex_tree.fill();
      continue;
    }

    // FIXME: better title for canvas
    std::string canvas_title{"event"};
    for (const auto& path : paths) {
      canvas_title.push_back(' ');
      canvas_title.insert(canvas_title.cend(), path.begin(), path.end());
    }

    plot::canvas canvas("event" + event_counter_string,
                        canvas_title + " | " + event_counter_string,
                        !options.verbose_output);
    if (options.draw_events) {
      box::io::draw_detector(canvas, extension.layer_count);
      box::io::draw_mc_tracks(canvas, mc::convert_events(mc_imported_events[event_counter]));
      for (const auto& hit : altered_event)
        canvas.add_point(type::reduce_to_r3(hit), 0.8, plot::color::BLACK);
    }

    analysis::full_event non_primary_track_points, non_secondary_track_points;
    auto tracks = find_tracks(altered_event,
                              options,
                              20.0L,
                              1UL,
                              canvas,
                              non_primary_track_points);
    auto secondary_tracks = find_tracks(non_primary_track_points,
                                        options,
                                        100.0L,
                                        0UL,
                                        canvas,
                                        non_secondary_track_points);
    tracks.reserve(tracks.size() + secondary_tracks.size());
    tracks.insert(tracks.cend(),
                  std::make_move_iterator(secondary_tracks.cbegin()),
                  std::make_move_iterator(secondary_tracks.cend()));

    box::io::save_tracks(tracks, canvas, track_tree, options);
    box::io::print_tracking_summary(event, tracks);

    analysis::track_vector converged_tracks;
    converged_tracks.reserve(tracks.size());
    std::copy_if(std::make_move_iterator(tracks.begin()),
                 std::make_move_iterator(tracks.end()),
                 std::back_inserter(converged_tracks),
                 [](const auto& track) { return track.fit_converged(); });

    box::io::save_vertices(analysis::pairwise_fit_tracks(converged_tracks), canvas, vertex_tree, options);

    canvas.draw();
  }

  plot::save_all(save_path,
    filetype_tag,
    plot::value_tag{"TIMESTAMP", util::time::GetString("%c %Z")},
    project_tag,
    plot::value_tag{"EVENTS", std::to_string(import_size)},
    plot::value_tag{"EFFICIENCY", std::to_string(options.simulated_efficiency)},
    plot::value_tag{"NOISE", std::to_string(options.simulated_noise_rate * units::time) + " / " + units::time_string},
    box::geometry::value_tags(),
    box::io::data_paths_value_tags(paths, options.data_timing_offsets));



  const auto path_count = save_path.size();
  std::vector<std::string> prefixes;
  prefixes.reserve(path_count);
  for (std::size_t i{}; i < path_count; ++i)
      prefixes.push_back("SIM_" + std::to_string(i) + "_");

  if (options.merge_input) {
      reader::root::merge_save(save_path, paths, prefixes);
	  for (std::size_t i{}; i < path_count; ++i) {
          track_tree.add_friend(prefixes[i] + "box_run", save_path);
          vertex_tree.add_friend(prefixes[i] + "box_run", save_path);
          const std::string name = prefixes[i] + "box_run";
          digi_tree.AddFriend(name.c_str(), save_path.c_str());
	  }
  }

  track_tree.save(save_path);
  vertex_tree.save(save_path);

  auto file = TFile::Open(save_path.c_str(), "UPDATE");
  if (file && !file->IsZombie()) {
	  file->cd();
	  digi_tree.Write();
	  file->Close();
  }



}
//----------------------------------------------------------------------------------------------

//__Box Tracking Algorithm______________________________________________________________________
int box_tracking(int argc,
                 char* argv[]) {
  box::io::extension_parser extension{};
  const auto options = script::parse_command_line(argc, argv, extension);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file, options.default_time_error);
  box::geometry::import(extension);

  box::io::print_tracking_directories(options.data_directories);

  std::size_t path_counter{};
  const auto statistics_path_prefix = box::io::add_statistics_path(options);
  for (const auto& paths : reader::root::transpose_search_directories(options.data_directories, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);
    const auto statistics_save_path = statistics_path_prefix
                                    + path_counter_string
                                    + "."
                                    + options.statistics_file_extension;

    box::io::print_bar();
    box::io::print_tracking_paths(paths);

    track_event_bundle(paths,
      reader::root::import_full_event_mc_bundle(paths,
          options.data_timing_offsets,
          options.data_track_id_key,
          options.data_t_key,
          options.data_x_key,
          options.data_y_key,
          options.data_z_key,
          options.data_e_key,
          options.data_px_key,
          options.data_py_key,
          options.data_pz_key,
          options.data_detector_key),
      options,
      extension,
      statistics_save_path);

    std::cout << "\n\nSaved: " << statistics_save_path << "\n";
  }

  box::io::print_bar();
  geometry::close();
  plot::end();
  std::cout << "Completed.\n";
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Silent Box Tracking Algorithm_______________________________________________________________
int silent_box_tracking(int argc,
                        char* argv[]) {
  util::io::remove_buffer(std::cout, std::cerr, std::clog);
  return box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------



} /* namespace MATHUSLA */

//__Main Function: Box Tracker__________________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------
