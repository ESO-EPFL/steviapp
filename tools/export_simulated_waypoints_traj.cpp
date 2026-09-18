
#include "datablocks/project.h"
#include "datablocks/trajectory.h"

#include "geo/wgs84.h"

#include "testutils/eceftrajectoryinertialinstrumentssimulator.h"
#include "testutils/waypointssimulatedtrajectorysmoother.h"

#include <memory>

#include <QCommandLineParser>
#include <QTextStream>
#include <QFile>
#include <QFileInfo>
#include <QDir>

int main (int argc, char** argv) {

    QTextStream out(stdout);
    QTextStream err(stderr);

    QCommandLineParser parser;

    parser.setApplicationDescription("Export simulated trajectory from waypoints");
    QString out_proj_option_name = "out_project";
    parser.addPositionalArgument(out_proj_option_name,
                                 "output project file. Additional files will be written in the same directory.");

    QString waypoints_option_name = "waypoints";
    parser.addPositionalArgument(waypoints_option_name,
                                 "file containing the waypoints of the trajectory.");

    QString gpsOptionName = "gpsFreq";
    QString insOptionName = "insFreq";
    QString speedOptionName = "speed";
    QString smoothOptionName = "smoothTime";
    QCommandLineOption gpsFreqOption(QStringList{"g","gps_freq"}, "frequency of the gps sampling, in Hz.", gpsOptionName, "2");
    QCommandLineOption insFreqOption(QStringList{"i","ins_freq"}, "frequency of the ins sampling, in Hz.", insOptionName, "10");
    QCommandLineOption speedOption(QStringList{"s","speed"}, "frequency of the gps sampling, in Hz.", speedOptionName, "50");
    QCommandLineOption smoothTimeOption(QStringList{"t","smooth_time"}, "frequency of the ins sampling, in Hz.", smoothOptionName, "4");
    parser.addOptions({gpsFreqOption, insFreqOption, speedOption, smoothTimeOption});

    parser.addHelpOption();

    QStringList args;
    args.reserve(argc);

    for (int i = 0; i < argc; i++) {
        args << argv[i];
    }

    parser.parse(args);

    QStringList posArgs = parser.positionalArguments();

    QString out_proj_path = posArgs[0];
    QString in_waypoints_path = posArgs[1];

    QFile in_waypoints(in_waypoints_path);

    if (!in_waypoints.open(QFile::ReadOnly)) {
        err << "Error opening waypoints file" << Qt::endl;
        return 1;
    }

    int freqGps = 2; // gps sampling freq Hz
    int freqIns = 10; // ins sampling freq Hz

    double speed = 10; // 50 m/s
    double smoothTime = 8; // 0.5 s

    QString gpsArg = parser.value(gpsFreqOption);

    if (!gpsArg.isEmpty()) {
        bool ok = true;
        freqGps = gpsArg.toInt(&ok);

        if (!ok or freqGps <= 0) {
            err << "Error parsing gps freq argument" << Qt::endl;
            return 1;
        }
    }

    QString insArg = parser.value(insFreqOption);

    if (!insArg.isEmpty()) {
        bool ok = true;
        freqIns = insArg.toInt(&ok);

        if (!ok or freqIns <= 0) {
            err << "Error parsing ins freq argument" << Qt::endl;
            return 1;
        }
    }

    QString speedArg = parser.value(speedOption);

    if (!speedArg.isEmpty()) {
        bool ok = true;
        speed = speedArg.toDouble(&ok);

        if (!ok or speed <= 0) {
            err << "Error parsing speed argument" << Qt::endl;
            return 1;
        }
    }

    QString smoothTimeArg = parser.value(smoothTimeOption);

    if (!smoothTimeArg.isEmpty()) {
        bool ok = true;
        smoothTime = smoothTimeArg.toDouble(&ok);

        if (!ok or smoothTime <= 0) {
            err << "Error parsing smooth time argument" << Qt::endl;
            return 1;
        }
    }

    std::vector<std::array<double,3>> wayPointsWGS84;

    while (!in_waypoints.atEnd()) {
        QByteArray lineData = in_waypoints.readLine();
        QString line = QString::fromLocal8Bit(lineData);
        if (line.isEmpty() or line.startsWith("#")) {
            continue;
        }
        QStringList splitted = line.split(QRegExp("\\s"), Qt::SkipEmptyParts);

        if (splitted.size() != 3) {
            err << "Error while reading waypoints file" << Qt::endl;
            return 1;
        }

        double lat;
        double lon;
        double alt;

        bool ok = true;

        lat = splitted[0].toDouble(&ok);

        if (!ok) {
            err << "Error while reading waypoints file" << Qt::endl;
            return 1;
        }

        lon = splitted[1].toDouble(&ok);

        if (!ok) {
            err << "Error while reading waypoints file" << Qt::endl;
            return 1;
        }

        alt = splitted[2].toDouble(&ok);

        if (!ok) {
            err << "Error while reading waypoints file" << Qt::endl;
            return 1;
        }

        wayPointsWGS84.push_back({lat, lon, alt});
    }

    int tickGPS = freqIns/freqGps;

    double dtPos = 1./freqGps;
    double dtAcc = 1./freqIns;
    double dt = std::min(dtAcc, dtPos);

    if (dt <= 0) {
        err << "Error in time configuration" << Qt::endl;
        return 1;
    }

    std::vector<std::array<double,3>> waypointsECEF(wayPointsWGS84.size());

    Eigen::Vector3d meanPos = Eigen::Vector3d::Zero();

    for (size_t i = 0; i < wayPointsWGS84.size(); i++) {
        waypointsECEF[i] = StereoVisionApp::Geo::WGS84Ellipsoid::LatLonHeight2ECEF(wayPointsWGS84[i]);
        for (int j = 0; j < 3; j++) {
            meanPos[j] += waypointsECEF[i][j];
        }
    }

    meanPos /= waypointsECEF.size();

    out << "wgs84 landmarks: \n";
    for (std::array<double,3> const& lm : wayPointsWGS84) {
        out << "\t" << lm[0] << " " << lm[1] << " " << lm[2] << "\n";
    }
    out << "ecef landmarks: \n";
    for (std::array<double,3> const& lm : waypointsECEF) {
        out << "\t" << lm[0] << " " << lm[1] << " " << lm[2] << "\n";
    }
    out << "Mean pos ecef: " << meanPos.x() << " " << meanPos.y() << " " << meanPos.z() << Qt::endl;

    StereoVisionApp::Simulation::WayPointsEcefTrajectoryFunctor* trajFunctor =
        new StereoVisionApp::Simulation::WayPointsEcefTrajectoryFunctor(waypointsECEF,speed,smoothTime);

    StereoVisionApp::Simulation::EcefTrajectoryInertialInstrumentsSimulator staticTrajSimulator(trajFunctor);

    auto wayPointsInternal = trajFunctor->orientationWaypoints();
    double t0 = wayPointsInternal.times.front();
    double tf = wayPointsInternal.times.back();

    int nSamples = std::ceil((tf-t0)/dt);
    int nSamplesGps = std::ceil((tf-t0)/dtPos);

    double currentT = t0;
    int tickCount = 0;

    std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock> posCache;
    std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock> orientCache;
    StereoVisionApp::Trajectory::RawGpsData gpsCache;
    std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock> gyroCache;
    std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock> accCache;

    posCache.reserve(nSamplesGps);
    orientCache.reserve(nSamplesGps);
    gpsCache.position = std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock>();
    gpsCache.position->reserve(nSamplesGps);
    gpsCache.velocities = std::vector<StereoVisionApp::Trajectory::TimeCartesianBlock>();
    gpsCache.velocities->reserve(nSamplesGps);

    gyroCache.reserve(nSamples);
    accCache.reserve(nSamples);

    const StereoVisionApp::Geo::TopocentricConvention topocentricConvention =
        (StereoVisionApp::Geo::TopocentricConvention)(StereoVisionApp::Simulation::EcefTrajectoryInertialInstrumentsSimulator::LocalFrameDefinitionUsed);

    while (currentT < tf) {

        if (tickCount % tickGPS == 0) {
            StereoVision::Geometry::RigidBodyTransform<double> body2ecef = staticTrajSimulator.body2ecef(currentT);

            posCache.push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=body2ecef.t});

            Eigen::Matrix3d topocentric2ecef = StereoVisionApp::Geo::localFrame2ECEFFromECEF(body2ecef.t, topocentricConvention);
            Eigen::Vector3d orient = StereoVision::Geometry::inverseRodriguezFormula<double>(topocentric2ecef.transpose()*
                                                                                             StereoVision::Geometry::rodriguezFormula<double>(body2ecef.r)); //body2topocentric

            orientCache.push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=orient});

            Eigen::Vector3d gps = staticTrajSimulator.gps(currentT);
            Eigen::Vector3d gpsVelocity = staticTrajSimulator.gpsVelocity(currentT);

            gpsCache.position->push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=gps});
            gpsCache.velocities->push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=gpsVelocity});
        }
        Eigen::Vector3d gyro = staticTrajSimulator.gyro(currentT);
        Eigen::Vector3d acc = staticTrajSimulator.acc(currentT);

        gyroCache.push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=gyro});
        accCache.push_back(StereoVisionApp::Trajectory::TimeCartesianBlock{.time=currentT,.val=acc});

        currentT += dt;
        tickCount++;
    }

    QFileInfo projFileInfo(out_proj_path);
    QDir outDir = projFileInfo.dir();

    QString trajFilePath = outDir.absoluteFilePath("trajectory.csv");
    QString insFilePath = outDir.absoluteFilePath("ins.csv");

    QFile trajFile(trajFilePath);
    trajFile.open(QFile::WriteOnly);
    if (!trajFile.isOpen()) {
        err << "error opening output trajectory file" << Qt::endl;
        return 1;
    }

    QTextStream outTraj(&trajFile);

    outTraj << "#pos is in ecef (epsg 4978), orientation w.r.t. local frame ENU and expressed as axis angle, gps is epsg 4978\n";
    outTraj << "#time,pos_x,pos_y,pos_z,orient_x,orient_y,orient_z,gps_x,gps_y,gps_z,gps_vel_x,gps_vel_y,gps_vel_z\n";
    outTraj.setRealNumberPrecision(4);
    outTraj.setRealNumberNotation(QTextStream::RealNumberNotation::FixedNotation);
    for (size_t i = 0; i < posCache.size(); i++) {

        if (posCache[i].time != orientCache[i].time) {
            err << "timing error between positions and orientations, aborting!" << Qt::endl;
            return 1;
        }

        if (posCache[i].time != gpsCache.position->at(i).time) {
            err << "timing error between positions and gps, aborting!" << Qt::endl;
            return 1;
        }

        if (posCache[i].time != gpsCache.velocities->at(i).time) {
            err << "timing error between positions and gps velocities, aborting!" << Qt::endl;
            return 1;
        }

        outTraj << posCache[i].time << ", ";

        outTraj << posCache[i].val[0] << ", ";
        outTraj << posCache[i].val[1] << ", ";
        outTraj << posCache[i].val[2] << ", ";

        outTraj << orientCache[i].val[0] << ", ";
        outTraj << orientCache[i].val[1] << ", ";
        outTraj << orientCache[i].val[2] << ", ";

        outTraj << gpsCache.position->at(i).val[0] << ", ";
        outTraj << gpsCache.position->at(i).val[1] << ", ";
        outTraj << gpsCache.position->at(i).val[2] << ", ";

        outTraj << gpsCache.velocities->at(i).val[0] << ", ";
        outTraj << gpsCache.velocities->at(i).val[1] << ", ";
        outTraj << gpsCache.velocities->at(i).val[2] << "\n";
    }

    trajFile.close();


    QFile insFile(insFilePath);
    insFile.open(QFile::WriteOnly);
    if (!insFile.isOpen()) {
        err << "error opening output inertial data file" << Qt::endl;
        return 1;
    }

    QTextStream outIns(&insFile);

    outIns << "#specific force in local frame, gyro uses axis angle representation\n";
    outIns << "#time,acc_x,acc_y,acc_z,gyro_x,gyro_y,gyro_z\n";
    for (size_t i = 0; i < accCache.size(); i++) {

        if (accCache[i].time != gyroCache[i].time) {
            err << "timing error between accelerometer and gyro, aborting!" << Qt::endl;
            return 1;
        }

        outIns << accCache[i].time << ", ";

        outIns << accCache[i].val[0] << ", ";
        outIns << accCache[i].val[1] << ", ";
        outIns << accCache[i].val[2] << ", ";

        outIns << gyroCache[i].val[0] << ", ";
        outIns << gyroCache[i].val[1] << ", ";
        outIns << gyroCache[i].val[2] << "\n";

    }

    insFile.close();

    StereoVisionApp::ProjectFactory& pF = StereoVisionApp::ProjectFactory::defaultProjectFactory();
    pF.addType(new StereoVisionApp::TrajectoryFactory(&pF));

    std::unique_ptr<StereoVisionApp::Project> pPtr(pF.createProject());

    if (pPtr == nullptr) {
        err << "error creating project!" << Qt::endl;
        return 1;
    }

    StereoVisionApp::Project& project = *pPtr;
    const char* ecefCRS = "EPSG:4978";

    project.setDefaultProjectCRS(ecefCRS);

    qint64 trajectoryId = project.createDataBlock(StereoVisionApp::Trajectory::staticMetaObject.className());

    StereoVisionApp::Trajectory* traj = project.getDataBlock<StereoVisionApp::Trajectory>(trajectoryId);

    if (traj == nullptr) {
        err << "error creating trajectory block" << Qt::endl;
        return 1;
    }

    //position
    traj->setPositionFile(trajFilePath);
    traj->setPositionColumn(StereoVisionApp::Trajectory::Axis::T,0);
    traj->setPositionColumn(StereoVisionApp::Trajectory::Axis::X,1);
    traj->setPositionColumn(StereoVisionApp::Trajectory::Axis::Y,2);
    traj->setPositionColumn(StereoVisionApp::Trajectory::Axis::Z,3);

    traj->setPositionEpsg(ecefCRS);

    traj->setPositionTimeDelta(0);
    traj->setPositionTimeScale(1);

    //orientation
    traj->setOrientationFile(trajFilePath);
    traj->setOrientationColumn(StereoVisionApp::Trajectory::Axis::T,0);
    traj->setOrientationColumn(StereoVisionApp::Trajectory::Axis::X,4);
    traj->setOrientationColumn(StereoVisionApp::Trajectory::Axis::Y,5);
    traj->setOrientationColumn(StereoVisionApp::Trajectory::Axis::Z,6);

    traj->setOrientationSign(StereoVisionApp::Trajectory::Axis::X,1);
    traj->setOrientationSign(StereoVisionApp::Trajectory::Axis::Y,1);
    traj->setOrientationSign(StereoVisionApp::Trajectory::Axis::Z,1);

    traj->setOrientationTimeDelta(0);
    traj->setOrientationTimeScale(1);

    traj->setOrientationAngleRepresentation(StereoVisionApp::Trajectory::AxisAngle);
    traj->setOrientationAngleUnits(StereoVisionApp::Trajectory::Radians);
    traj->setOrientationTopocentricConvention(topocentricConvention);

    //accelerometer
    traj->setAccelerometerFile(insFilePath);
    traj->setAccelerometerColumn(StereoVisionApp::Trajectory::Axis::T,0);
    traj->setAccelerometerColumn(StereoVisionApp::Trajectory::Axis::X,1);
    traj->setAccelerometerColumn(StereoVisionApp::Trajectory::Axis::Y,2);
    traj->setAccelerometerColumn(StereoVisionApp::Trajectory::Axis::Z,3);

    traj->setAccelerometerTimeDelta(0);
    traj->setAccelerometerTimeScale(1);

    //gyroscope
    traj->setGyroFile(insFilePath);
    traj->setGyroColumn(StereoVisionApp::Trajectory::Axis::T,0);
    traj->setGyroColumn(StereoVisionApp::Trajectory::Axis::X,4);
    traj->setGyroColumn(StereoVisionApp::Trajectory::Axis::Y,5);
    traj->setGyroColumn(StereoVisionApp::Trajectory::Axis::Z,6);

    traj->setGyroAngleRepresentation(StereoVisionApp::Trajectory::AxisAngle);
    traj->setGyroAngleUnits(StereoVisionApp::Trajectory::Radians);

    traj->setGyroTimeDelta(0);
    traj->setGyroTimeScale(1);

    traj->setGyroSign(StereoVisionApp::Trajectory::Axis::X,1);
    traj->setGyroSign(StereoVisionApp::Trajectory::Axis::Y,1);
    traj->setGyroSign(StereoVisionApp::Trajectory::Axis::Z,1);

    //gps
    traj->setGpsFile(trajFilePath);
    traj->setGpsTopocentricConvention(topocentricConvention);
    traj->setGpsEpsg(ecefCRS);

    traj->setGpsTimeDelta(0);
    traj->setGpsTimeScale(1);

    traj->setGpsPosColumn(StereoVisionApp::Trajectory::Axis::T,0);
    traj->setGpsPosColumn(StereoVisionApp::Trajectory::Axis::X,7);
    traj->setGpsPosColumn(StereoVisionApp::Trajectory::Axis::Y,8);
    traj->setGpsPosColumn(StereoVisionApp::Trajectory::Axis::Z,9);

    traj->setGpsSpeedColumn(StereoVisionApp::Trajectory::Axis::X,10);
    traj->setGpsSpeedColumn(StereoVisionApp::Trajectory::Axis::Y,11);
    traj->setGpsSpeedColumn(StereoVisionApp::Trajectory::Axis::Z,12);

    project.save(out_proj_path);

    return 0;
}
