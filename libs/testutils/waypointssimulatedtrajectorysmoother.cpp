#include "./waypointssimulatedtrajectorysmoother.h"

#include "geo/wgs84.h"

namespace StereoVisionApp {
namespace Simulation {


WayPointsEcefTrajectoryFunctor WayPointsEcefTrajectoryFunctor::fromWGS84WayPoints(std::vector<std::array<double,3>> const& wayPointsWgs84,
                                                                                  double speed,
                                                                                  double smoothTime) {
    std::vector<std::array<double,3>> wayPointsEcef(wayPointsWgs84.size());

    for (size_t i = 0; i < wayPointsWgs84.size(); i++) {
        wayPointsEcef[i] = Geo::WGS84Ellipsoid::LatLonHeight2ECEF(wayPointsWgs84[i]);
    }

    return WayPointsEcefTrajectoryFunctor(wayPointsEcef, speed, smoothTime);
}

WayPointsEcefTrajectoryFunctor::WayPointsEcefTrajectoryFunctor(std::vector<std::array<double,3>> const& wayPointsEcef, double speed, double smoothTime) :
    _smoothTime(smoothTime)
{

    std::vector<double> times = WayPointsSimulatedTrajectorySmoother::estimateTimesFromWaypoints(wayPointsEcef, speed);

    std::array<double,3> mean_pos{0,0,0};

    for (size_t i = 0; i < wayPointsEcef.size(); i++) {

        for (int j = 0; j < 3; j++) {
            mean_pos[j] += wayPointsEcef[i][j];
        }
    }
    for (int j = 0; j < 3; j++) {
        mean_pos[j] /= wayPointsEcef.size();
    }

    std::array<double,3> gravity = Geo::WGS84Ellipsoid::gravityEcefModel(mean_pos);

    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> orientations =
        WayPointsSimulatedTrajectorySmoother::estimateAircraftOrientation(wayPointsEcef,times, gravity, smoothTime);

    _position_waypoints.times.resize(times.size());
    _position_waypoints.samples.resize(wayPointsEcef.size());

    for (size_t i = 0; i < wayPointsEcef.size(); i++) {
        _position_waypoints.times[i] = times[i];
        for (int j = 0; j < 3; j++) {
            _position_waypoints.samples[i][j] = wayPointsEcef[i][j];
        }
    }

    _orientation_waypoints.times.resize(orientations.times.size());
    _orientation_waypoints.samples.resize(orientations.samples.size());

    for (size_t i = 0; i < orientations.samples.size(); i++) {
        _orientation_waypoints.times[i] = orientations.times[i];
        for (int j = 0; j < 3; j++) {
            _orientation_waypoints.samples[i][j] = orientations.samples[i][j];
        }
    }

}

WayPointsEcefTrajectoryFunctor::~WayPointsEcefTrajectoryFunctor() {

}


StereoVision::Geometry::RigidBodyTransform<ceres::Jet<ceres::Jet<double,1>,1>> WayPointsEcefTrajectoryFunctor::trajectory(ceres::Jet<ceres::Jet<double,1>,1> const& t_ext) const {

    D2Jet t = t_ext;

    std::array<D2Jet,3> pos = WayPointsSimulatedTrajectorySmoother::sampleTrajectory(t, _position_waypoints.samples, _position_waypoints.times, _smoothTime);
    std::array<D2Jet,3> rot = WayPointsSimulatedTrajectorySmoother::sampleTrajectory(t, _orientation_waypoints.samples, _orientation_waypoints.times, _smoothTime);

    return StereoVision::Geometry::RigidBodyTransform<ceres::Jet<ceres::Jet<double,1>,1>>(
        Eigen::Matrix<ceres::Jet<ceres::Jet<double,1>,1>,3,1>(rot[0],rot[1],rot[2]),
        Eigen::Matrix<ceres::Jet<ceres::Jet<double,1>,1>,3,1>(pos[0],pos[1],pos[2]));

}

WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> WayPointsEcefTrajectoryFunctor::positionWaypoints() const {
    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> ret;
    ret.times.reserve(_position_waypoints.times.size());
    ret.samples.reserve(_position_waypoints.samples.size());

    for (int i = 0; i < _position_waypoints.times.size(); i++) {
        ret.times.push_back(double(_position_waypoints.times[i]));
    }

    for (int i = 0; i < _position_waypoints.samples.size(); i++) {
        std::array<double,3> sample;
        for (int j = 0; j < 3; j++) {
            sample[j] = double(_position_waypoints.samples[i][j]);
        }
        ret.samples.push_back(sample);
    }

    return ret;
}
WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> WayPointsEcefTrajectoryFunctor::orientationWaypoints() const {
    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> ret;
    ret.times.reserve(_orientation_waypoints.times.size());
    ret.samples.reserve(_orientation_waypoints.samples.size());

    for (int i = 0; i < _orientation_waypoints.times.size(); i++) {
        ret.times.push_back(double(_orientation_waypoints.times[i]));
    }

    for (int i = 0; i < _orientation_waypoints.samples.size(); i++) {
        std::array<double,3> sample;
        for (int j = 0; j < 3; j++) {
            sample[j] = double(_orientation_waypoints.samples[i][j]);
        }
        ret.samples.push_back(sample);
    }

    return ret;
}

} // namespace Simulation
} // namespace StereoVisionApp
