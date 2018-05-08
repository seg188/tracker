#ifndef TRACKER__READER_HH
#define TRACKER__READER_HH
#pragma once

#include <unordered_map>

#include "analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

using namespace type;

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__ROOT Directory Search_______________________________________________________________________
std::vector<std::string> search_directory(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Detector Map___________________________________________________________________________
using detector_map = std::unordered_map<type::integer, std::string>;
//----------------------------------------------------------------------------------------------

//__ROOT Detector Map Import____________________________________________________________________
detector_map import_detector_map(const std::string& path);
//----------------------------------------------------------------------------------------------

//__ROOT Event Import___________________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& x_key,
                                     const std::string& y_key,
                                     const std::string& z_key);
analysis::event_vector import_events(const std::string& path,
                                     const std::string& time_key,
                                     const std::string& detector_key,
                                     const detector_map& map);
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

//__Tracking Script Options_____________________________________________________________________
struct tracking_options {
  std::string geometry_file;
  std::string geometry_map_file;
  std::string root_directory;
  std::string root_time_key;
  std::string root_x_key;
  std::string root_y_key;
  std::string root_z_key;
  std::string root_detector_key;
  r4_point    collapse_size;
  real        layer_depth;
  real        line_width;
  integer     seed_size;
};
const static tracking_options default_options{
  "", "", "", "Time", "X", "Y", "Z", "Detector", {0, 0, 0, 0}, 500, 1, 3};
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path);
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__READER_HH */
