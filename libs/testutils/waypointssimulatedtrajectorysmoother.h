#ifndef WAYPOINTSSIMULATEDTRAJECTORYSMOOTHER_H
#define WAYPOINTSSIMULATEDTRAJECTORYSMOOTHER_H

#include <array>
#include <vector>

#include <StereoVision/geometry/rotations.h>

#include "./eceftrajectoryinertialinstrumentssimulator.h"

namespace StereoVisionApp {
namespace Simulation {

class WayPointsSimulatedTrajectorySmoother {
public:

    template<typename T>
    struct TimesWaypointsSet {
        std::vector<T> times;
        std::vector<std::array<T,3>> samples;
    };

    template <typename T, typename WpContT>
    static std::vector<T> estimateTimesFromWaypoints(WpContT const& wayPoints, T speed) {
        std::vector<T> times(wayPoints.size());
        times[0] = 0;

        for (int i = 1; i < wayPoints.size(); i++) {
            T d = 0;
            for (int j = 0; j < 3; j++) {
                T tmp = wayPoints[i][j] - wayPoints[i-1][j];
                d += tmp*tmp;
            }
            d = sqrt(d);

            T t0(0);
            if (i > 0) {
                t0 = times[i-1];
            }
            times[i] = t0 + d/speed;
        }

        return times;
    }

    template <typename T, typename ContT, typename WpContT>
    static TimesWaypointsSet<T> estimateAircraftOrientation(WpContT const& wayPoints, ContT const& ts, std::array<T,3> const& gravity, T smoothingKernelRadius) {
        int nTrajWayPoints = wayPoints.size();

        if (nTrajWayPoints <= 0) {
            return TimesWaypointsSet<T>{.times={},.samples={}};
        }

        if (nTrajWayPoints == 1) {
            return TimesWaypointsSet<T>{.times={ts[0]},.samples={std::array<T,3>{T(0),T(0),T(0)}}};
        }

        Eigen::Matrix<T,3,1> g_vec;
        g_vec[0] = -gravity[0];
        g_vec[1] = -gravity[1];
        g_vec[2] = -gravity[2];
        Eigen::Matrix<T,3,1> down = g_vec.normalized();

        auto computeOrientationForSegment = [&down] (Eigen::Matrix<T,3,1> const& segmentDir) -> std::array<T,3> {
            T norm = segmentDir.norm();
            if (norm < 1e-6) {
                return std::array<T,3>{T(0),T(0),T(0)};
            }
            Eigen::Matrix<T,3,1> n_x = segmentDir/norm;
            Eigen::Matrix<T,3,1> n_y = down.cross(n_x);
            T normY = n_y.norm();
            if (normY < 1e-6) {
                Eigen::Matrix<T,3,1> x_vec(T(1),T(0),T(0));
                Eigen::Matrix<T,3,1> d = x_vec.cross(n_x);
                T normD = d.norm();
                if (normD < 1e-6) {
                    return std::array<T,3>{T(0),T(0),T(0)};
                }
                T theta = std::acos(normD);
                d /= normD;
                d *= theta;

                return std::array<T,3>{d[0],d[1],d[2]};
            }
            n_y /= normY;
            Eigen::Matrix<T,3,1> n_z = n_x.cross(n_y);

            Eigen::Matrix<T,3,3> R;
            R.col(0) = n_x;
            R.col(1) = n_y;
            R.col(2) = n_z;

            Eigen::Matrix<T,3,1> r = StereoVision::Geometry::inverseRodriguezFormula(R);
            return std::array<T,3>{r[0],r[1],r[2]};
        };

        if (nTrajWayPoints == 2) {
            Eigen::Matrix<T,3,1> segmentDir;
            for (int i = 0; i < 3; i++) {
                segmentDir[i] = wayPoints[1][i] - wayPoints[0][i];
            }
            return TimesWaypointsSet<T>{.times={ts[0]},.samples={computeOrientationForSegment(segmentDir)}};
        }

        auto computeOrientationBetweenSegment = [&g_vec, &computeOrientationForSegment, &smoothingKernelRadius] (Eigen::Matrix<T,3,1> const& segment1Speed,
                                                                                        Eigen::Matrix<T,3,1> const& segment2Speed) -> std::array<T,3> {
            T acc_dt = smoothingKernelRadius*T(2);
            Eigen::Matrix<T,3,1> acc = (segment2Speed - segment1Speed)/acc_dt;
            Eigen::Matrix<T,3,1> midSpeed = (segment2Speed + segment1Speed)/2;
            T acc_norm = acc.norm();
            if (acc_norm < 1e-6) {
                return computeOrientationForSegment(midSpeed);
            }
            T dir_norm = midSpeed.norm();
            if (dir_norm < 1e-6) {
                return computeOrientationForSegment(midSpeed);
            }

            Eigen::Matrix<T,3,1> modif_down = g_vec-acc;

            Eigen::Matrix<T,3,1> n_x = midSpeed/dir_norm;
            Eigen::Matrix<T,3,1> n_y = modif_down.cross(n_x);
            T normY = n_y.norm();
            if (normY < 1e-6) {
                Eigen::Matrix<T,3,1> x_vec(T(1),T(0),T(0));
                Eigen::Matrix<T,3,1> d = x_vec.cross(n_x);
                T normD = d.norm();
                if (normD < 1e-6) {
                    return std::array<T,3>{T(0),T(0),T(0)};
                }
                T theta = std::acos(normD);
                d /= normD;
                d *= theta;

                return std::array<T,3>{d[0],d[1],d[2]};
            }
            n_y /= normY;
            Eigen::Matrix<T,3,1> n_z = n_x.cross(n_y);

            Eigen::Matrix<T,3,3> R;
            R.col(0) = n_x;
            R.col(1) = n_y;
            R.col(2) = n_z;

            Eigen::Matrix<T,3,1> r = StereoVision::Geometry::inverseRodriguezFormula(R);
            return std::array<T,3>{r[0],r[1],r[2]};
        };

        if (nTrajWayPoints == 2) {
            Eigen::Matrix<T,3,1> segmentDir;
            for (int i = 0; i < 3; i++) {
                segmentDir[i] = wayPoints[1][i] - wayPoints[0][i];
            }
            return TimesWaypointsSet<T>{.times={ts[0]},.samples={computeOrientationForSegment(segmentDir)}};
        }

        int nOrientWayPoints = 3*nTrajWayPoints - 4;

        TimesWaypointsSet<T> ret;
        ret.times.resize(nOrientWayPoints);
        ret.samples.resize(nOrientWayPoints);

        int r_i = 0;

        for (int i = 0; i < nTrajWayPoints-1; i++) {
            Eigen::Matrix<T,3,1> segmentDir;
            for (int j = 0; j < 3; j++) {
                segmentDir[j] = wayPoints[i+1][j] - wayPoints[i][j];
            }

            T dt = ts[i+1] - ts[i];

            std::array<T,3> r_segment = computeOrientationForSegment(segmentDir);

            ret.samples[r_i] = r_segment;
            T t = (i == 0) ? ts[i] : std::min(ts[i]+smoothingKernelRadius, ts[i+1]-smoothingKernelRadius);
            ret.times[r_i] = t;
            r_i++;

            ret.samples[r_i] = r_segment;
            t = (i == nTrajWayPoints-2) ? ts[i+1] : std::max(ts[i]+smoothingKernelRadius, ts[i+1]-smoothingKernelRadius);
            ret.times[r_i] = t;
            r_i++;

            if (i == nTrajWayPoints-2) {
                continue;
            }

            T next_dt = ts[i+2] - ts[i+1];

            Eigen::Matrix<T,3,1> nextSegmentDir;
            for (int j = 0; j < 3; j++) {
                nextSegmentDir[j] = wayPoints[i+2][j] - wayPoints[i+1][j];
            }

            std::array<T,3> next_r_segment = computeOrientationBetweenSegment(segmentDir/dt,nextSegmentDir/next_dt);

            ret.samples[r_i] = next_r_segment;
            t = ts[i+1];
            ret.times[r_i] = t;
            r_i++;

        }

        return ret;
    }

    /*!
     * \brief sampleTrajectory sample a trajectory from a set of waypoints
     * \param t the time at which the trajectory should be sampled
     * \param wayPoints the waypoints
     * \param ts the sampling times corresponding to the waypoints
     * \param smoothingKernelRadius the smoothing kernel's radius, larger mean more smoothing.
     * \return trajectory at time t.
     */
    template <typename T, typename ContT, typename WpContT>
    static std::array<T,3> sampleTrajectory(T t, WpContT const& wayPoints, ContT const& ts, T smoothingKernelRadius) {

        if (wayPoints.size() != ts.size()) {
            return std::array<T,3>{T(0), T(0), T(0)};
        }

        if (wayPoints.size() <= 0) {
            return std::array<T,3>{T(0), T(0), T(0)};
        }

        if (wayPoints.size() == 1) {
            return std::array<T,3>{wayPoints[0][0], wayPoints[0][1], wayPoints[0][2]};
        }

        std::vector<int> idxs;

        T t0 = t-smoothingKernelRadius;
        T t1 = t+smoothingKernelRadius;

        for (int i = 0; i < wayPoints.size()-1; i++) {
            if (ts[i] <= t0 and ts[i+1] >= t0) {
                idxs.push_back(i);
            } else if (ts[i] >= t0 and ts[i] <= t1) {
                idxs.push_back(i);
            } else if (ts[i] >= t1) {
                idxs.push_back(i);
                break;
            }
        }
        if (wayPoints.size() >= 2) {
            if (ts[wayPoints.size()-2] <= t1) {
                idxs.push_back(wayPoints.size()-1);
            }
        }

        std::vector<T> s_xs(idxs.size());
        std::vector<T> s_ys(idxs.size());
        std::vector<T> s_zs(idxs.size());
        std::vector<T> s_ts(idxs.size());

        for (size_t i = 0; i < idxs.size(); i++) {
            s_xs[i] = T(wayPoints[idxs[i]][0]);
            s_ys[i] = T(wayPoints[idxs[i]][1]);
            s_zs[i] = T(wayPoints[idxs[i]][2]);
            s_ts[i] = T(ts[idxs[i]]);
        }

        std::array<T,3> ret;
        ret[0] = evaluateSmoothedCurve(t, s_xs, s_ts, smoothingKernelRadius);
        ret[1] = evaluateSmoothedCurve(t, s_ys, s_ts, smoothingKernelRadius);
        ret[2] = evaluateSmoothedCurve(t, s_zs, s_ts, smoothingKernelRadius);

        return ret;
    }

    /*!
     * \brief evaluateSmoothedCurve smoothly interpolate a curve given its nodes
     * \param t the time at which the curve needs to be interpolated
     * \param ys the value of the linearily interpolated curve at control times
     * \param ts the control times
     * \param smoothingKernelRadius the scaling radius of the interpolation kernel. Large radius mean more smooting.
     * \return the value at time t of the smoothed curve.
     *
     * The function use an internal smoothing kernel choosen for the simulation, the exact value of which
     * might change in future versions, so you should not rely on the behavior of this class stay consistent in the future
     * when using it to write tests.
     *
     * The function assume that all interval in the provided control points are within the radius of the kernel.
     * If not, the function will still iterate over all nodes. It is the responsability of the user to
     * provide only relevant nodes in the ys and ts arrays. Times are assumed to be sorted.
     */
    template <typename T, typename ContT>
    static T evaluateSmoothedCurve(T t, ContT const& ys, ContT const& ts, T smoothingKernelRadius) {
        //current implementation use a kernel k(t) of the form y(t) = -2*t^3 -3*t^2 + 1 for t in [-1,0] and y(t) = 2*t^3 -3*t^2 + 1 for t in [0,1]
        //alternative expression is f(x) = 2∙x^2∙abs(x) −3∙x^2 + 1

        T ret(0);

        int nNodes = ts.size();

        if (nNodes == 0) {
            return ret;
        }

        if (ys.size() != nNodes) {
            return ret;
        }

        if (smoothingKernelRadius <= 0) {
            return ret;
        }

        if (nNodes == 1) {
            ret = ys[0];
            return ret;
        }

        //current start of interval (in scaled time)
        T i_ts;
        //current end of interval (in scaled time)
        T i_tf;
        //current evaluation time (scaled)
        T i_t = t/smoothingKernelRadius;
        //invicate if we already skipped the midpoint of the kernel.
        T sgn(-1);

        auto integrateInterval = [&i_t, &i_ts, &i_tf, &sgn] (T a, T b) -> T {
            T ret(0);
            // int_ts^tf k(tau-t) * (a*tau + b) - [tau = taumt+t / taumt = tau - t] -> int_taumt_s^taumt_f k(taumt) * (a*(taumt+t) + b)
            T taumt_f = i_tf-i_t;
            T taumt_s = i_ts-i_t;
            // int_taumt_s^taumt_f sgn*2*(taumt)^3 * (a*(taumt+t) + b)
            ret += sgn*taumt_f*taumt_f*taumt_f*taumt_f*(T(2./5.)*a*taumt_f + T(2./4.)*(a*i_t+b));
            ret -= sgn*taumt_s*taumt_s*taumt_s*taumt_s*(T(2./5.)*a*taumt_s + T(2./4.)*(a*i_t+b));
            // int_taumt_s^taumt_f -3*(taumt)^2 * (a*(taumt+t) + b)
            ret += -taumt_f*taumt_f*taumt_f*(T(3./4.)*a*taumt_f + T(3./3.)*(a*i_t+b));
            ret -= -taumt_s*taumt_s*taumt_s*(T(3./4.)*a*taumt_s + T(3./3.)*(a*i_t+b));
            // int_taumt_s^taumt_f 1 * (a*(taumt+t) + b)
            ret += taumt_f*(T(1./2.)*a*taumt_f + a*i_t+b);
            ret -= taumt_s*(T(1./2.)*a*taumt_s + a*i_t+b);
            return ret;
        };

        for (int i = 0; i < nNodes-1; i++) {
            T t0 = ts[i];
            T t1 = ts[i+1];
            T y0 = ys[i];
            T y1 = ys[i+1];

            T dt = t1 - t0;

            T a = y1/dt-y0/dt;
            T b = y0*t1/dt-y1*t0/dt;

            a *= smoothingKernelRadius;

            if (i == 0 and t1 > t-smoothingKernelRadius) {
                t0 = t-smoothingKernelRadius;
            }

            if (t > t0) { //middle of kernel is within interval

                i_ts = (t0 > t-smoothingKernelRadius) ? t0 : t-smoothingKernelRadius;
                if (i == 0) {
                    i_ts = t-smoothingKernelRadius;
                }
                i_tf = (t1 < t) ? t1 : t;
                i_ts /= smoothingKernelRadius;
                i_tf /= smoothingKernelRadius;

                dt = i_tf-i_ts;

                sgn = T(-1);

                if (dt > 0) {
                    ret += integrateInterval(a, b);
                }

            }

            if (i == nNodes-2 and t1 < t+smoothingKernelRadius) {
                t1 = t+smoothingKernelRadius;
            }

            if (t1 > t) {

                i_ts = (t0 > t) ? t0 : t;
                i_tf = (t1 < t+smoothingKernelRadius) ? t1 : t+smoothingKernelRadius;
                if (i == nNodes-2) {
                    i_tf = t+smoothingKernelRadius;
                }
                i_ts /= smoothingKernelRadius;
                i_tf /= smoothingKernelRadius;

                dt = i_tf-i_ts;

                sgn = T(1);

                if (dt > 0) {
                    ret += integrateInterval(a, b);
                }
            }

        }

        return ret;
    }



};

class WayPointsEcefTrajectoryFunctor : public InertialTrajectoryFunctor {
public :

    using D2Jet = InertialTrajectoryFunctor::D2_Jet;

    static WayPointsEcefTrajectoryFunctor fromWGS84WayPoints(std::vector<std::array<double,3>> const& wayPointsWgs84,
                                                             double speed = 10,
                                                             double smoothTime = 3);

    WayPointsEcefTrajectoryFunctor(std::vector<std::array<double,3>> const& wayPointsEcef,
                                   double speed = 10,
                                   double smoothTime = 3);
    virtual ~WayPointsEcefTrajectoryFunctor();


    virtual StereoVision::Geometry::RigidBodyTransform<ceres::Jet<ceres::Jet<double,1>,1>> trajectory(ceres::Jet<ceres::Jet<double,1>,1> const& t) const override;

    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> positionWaypoints() const;
    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<double> orientationWaypoints() const;

protected:

    D2Jet _smoothTime;
    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<D2Jet> _position_waypoints;
    WayPointsSimulatedTrajectorySmoother::TimesWaypointsSet<D2Jet> _orientation_waypoints;
};

} // namespace Simulation
} // namespace StereoVisionApp

#endif // WAYPOINTSSIMULATEDTRAJECTORYSMOOTHER_H
