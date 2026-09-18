#ifndef ECEFTRAJECTORYINERTIALINSTRUMENTSSIMULATOR_H
#define ECEFTRAJECTORYINERTIALINSTRUMENTSSIMULATOR_H

#include <ceres/jet.h> //used for automatic differentiation

#include <StereoVision/geometry/rotations.h>

namespace StereoVisionApp {
namespace Simulation {

class InertialTrajectoryFunctor {
public :
    /*!
     * \brief The D2_Jet class allow to use the ceres jet class for second order derivative computations a bit more conviniently with additional conversions defined.
     */
    class D2_Jet : public ceres::Jet<ceres::Jet<double,1>,1> {
    public:
        inline D2_Jet() : ceres::Jet<ceres::Jet<double,1>,1>() {
            a.a = 0;
            a.v[0] = 0;
            v[0].a = 0;
            v[0].v[0] = 0;
        }
        inline D2_Jet(ceres::Jet<double,1> const& d1) : ceres::Jet<ceres::Jet<double,1>,1>(d1) {
            a.a = d1.a;
            a.v[0] = d1.v[0];
            v[0].a = d1.v[0];
            v[0].v[0] = 0;
        }
        inline D2_Jet(ceres::Jet<ceres::Jet<double,1>,1> const& d2) : ceres::Jet<ceres::Jet<double,1>,1>(d2) {
            a.a = d2.a.a;
            a.v[0] = d2.a.v[0];
            v[0].a = d2.v[0].a;
            v[0].v[0] = d2.v[0].v[0];
        }
        inline D2_Jet(D2_Jet const& other) : ceres::Jet<ceres::Jet<double,1>,1>(other) {
            a.a = other.a.a;
            a.v[0] = other.a.v[0];
            v[0].a = other.v[0].a;
            v[0].v[0] = other.v[0].v[0];
        }
        inline D2_Jet(double scalar) : ceres::Jet<ceres::Jet<double,1>,1>() {
            a.a = scalar;
            a.v[0] = 0;
            v[0].a = 0;
            v[0].v[0] = 0;
        }
        inline D2_Jet(int scalar) : ceres::Jet<ceres::Jet<double,1>,1>() {
            a.a = scalar;
            a.v[0] = 0;
            v[0].a = 0;
            v[0].v[0] = 0;
        }

        inline D2_Jet& operator=(D2_Jet const& other) {
            a.a = other.a.a;
            a.v[0] = other.a.v[0];
            v[0].a = other.v[0].a;
            v[0].v[0] = other.v[0].v[0];
            return *this;
        }

        inline D2_Jet& operator=(ceres::Jet<double,1> const& other) {
            a.a = other.a;
            a.v[0] = other.v[0];
            v[0].a = other.v[0];
            v[0].v[0] = 0;
            return *this;
        }

        inline D2_Jet& operator=(double const& constant) {
            a.a = constant;
            a.v[0] = 0;
            v[0].a =  0;
            v[0].v[0] = 0;
            return *this;
        }

        inline bool operator< (ceres::Jet<ceres::Jet<double,1>,1> const& d) const {
            return a.a < d.a.a;
        }
        inline bool operator<= (ceres::Jet<ceres::Jet<double,1>,1> const& d) const {
            return a.a <= d.a.a;
        }
        inline bool operator> (ceres::Jet<ceres::Jet<double,1>,1> const& d) const {
            return a.a > d.a.a;
        }
        inline bool operator>= (ceres::Jet<ceres::Jet<double,1>,1> const& d) const {
            return a.a >= d.a.a;
        }

        inline bool operator< (double d) const {
            return a.a < d;
        }
        inline bool operator<= (double d) const {
            return a.a <= d;
        }
        inline bool operator> (double d) const {
            return a.a > d;
        }
        inline bool operator>= (double d) const {
            return a.a >= d;
        }

        inline bool operator< (int d) const {
            return a.a < d;
        }
        inline bool operator<= (int d) const {
            return a.a <= d;
        }
        inline bool operator> (int d) const {
            return a.a > d;
        }
        inline bool operator>= (int d) const {
            return a.a >= d;
        }

        inline operator double() const {
            return a.a;
        }
    };

    virtual ~InertialTrajectoryFunctor();

    /*!
     * \brief trajectory compute the trajectory, nested ceres jets allow to compute derivatives up to the second order automatically
     * \param t the time
     * \return the pose (body2inertial) at time t
     *
     * The trajectory is assumed to be given in the ecef frame.
     */
    virtual StereoVision::Geometry::RigidBodyTransform<ceres::Jet<ceres::Jet<double,1>,1>> trajectory(ceres::Jet<ceres::Jet<double,1>,1> const& t) const = 0;
};

/*!
 * \brief The EcefTrajectoryInertialInstrumentsSimulator class, from a trajectory function, give noiseless measurements for instruments used in navigation systems.
 *
 * This class use caching, so it is more efficient to collect all measurements for a given time, than all the measure for a single instrument, and then the next instrument, ect.
 */
class EcefTrajectoryInertialInstrumentsSimulator
{
public:

    static const int LocalFrameDefinitionUsed;

    struct Measurement {
        StereoVision::Geometry::RigidBodyTransform<double> body2ecef;
        Eigen::Vector3d gps; //gps measurement
        Eigen::Vector3d gpsVelocity; //gps velocity measurement
        Eigen::Vector3d gyro; //gyro measurement
        Eigen::Vector3d acc; //accelerometer measurement
    };

    /*!
     * \brief EcefTrajectoryInertialInstrumentsSimulator build a simulator with a given trajectory
     * \param trajectory the functor representing the trajectory. The simulator will take ownership of it
     */
    explicit EcefTrajectoryInertialInstrumentsSimulator(InertialTrajectoryFunctor* trajectory,
                                                        StereoVision::Geometry::RigidBodyTransform<double> const& gps2body =
                                                        StereoVision::Geometry::RigidBodyTransform<double>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()),
                                                        StereoVision::Geometry::RigidBodyTransform<double> const& ins2body =
                                                        StereoVision::Geometry::RigidBodyTransform<double>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()));
    ~EcefTrajectoryInertialInstrumentsSimulator();

    inline StereoVision::Geometry::RigidBodyTransform<double> body2ecef(double t) {
        if (std::isfinite(t) and t != _cachedMeasurementTime) {
            cacheMeasurement(t);
        }
        return _cachedMeasurement.body2ecef;
    }
    inline Eigen::Vector3d gps(double t) {
        if (std::isfinite(t) and t != _cachedMeasurementTime) {
            cacheMeasurement(t);
        }
        return _cachedMeasurement.gps;
    }
    inline Eigen::Vector3d gpsVelocity(double t) {
        if (std::isfinite(t) and t != _cachedMeasurementTime) {
            cacheMeasurement(t);
        }
        return _cachedMeasurement.gpsVelocity;
    }
    inline Eigen::Vector3d gyro(double t) {
        if (std::isfinite(t) and t != _cachedMeasurementTime) {
            cacheMeasurement(t);
        }
        return _cachedMeasurement.gyro;
    }
    inline Eigen::Vector3d acc(double t) {
        if (std::isfinite(t) and t != _cachedMeasurementTime) {
            cacheMeasurement(t);
        }
        return _cachedMeasurement.acc;
    }

protected:

    void cacheMeasurement(double t);

    double _cachedMeasurementTime;
    Measurement _cachedMeasurement;
    InertialTrajectoryFunctor* _functor;
    StereoVision::Geometry::RigidBodyTransform<double> _gps2body;
    StereoVision::Geometry::RigidBodyTransform<double> _ins2body;
};

} // namespace Simulation
} // namespace StereoVisionApp

#endif // ECEFTRAJECTORYINERTIALINSTRUMENTSSIMULATOR_H
