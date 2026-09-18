add_executable(export_simulated_waypoints_traj export_simulated_waypoints_traj.cpp)
target_link_libraries(export_simulated_waypoints_traj StereoVision::stevi ${TESTUTILS_LIB} ${UTILS_LIB} Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Test ${DATA_LIB})
