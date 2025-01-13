/*------------------------------------------------------------------------------
* pntpos_LEO.c : standard positioning
*
*          Copyright (C) 2007-2015 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2025/01/11 1.6  publish velocity from LEO Satellite
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

// add by weisong
#include <algorithm>
// google eigen
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>

#include <ros/ros.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PointStamped.h>
// #include <novatel_msgs/INSPVAX.h>
#include <novatel_oem7_msgs/INSPVAX.h> // novatel_msgs/INSPVAX
// #include <novatel_msgs/BESTPOS.h>
#include <novatel_oem7_msgs/BESTPOS.h> // novatel_msgs/INSPVAX

#include "../../include/gnss_tools.h"
#include <nlosExclusion/GNSS_Raw_Array.h>
#include <nlosExclusion/GNSS_Raw.h>

FILE* gnss_ublox_wls = fopen("gnss_ublox_wls.csv", "w+");


static const char rcsid[]="$Id:$";

/* constants -----------------------------------------------------------------*/
#define NFREQ       3           /* number of carrier frequencies */
#define SQR(x)      ((x)*(x))

#define NX          (4+3)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

ros::Publisher pub_pntpos_odometry;
ros::Publisher pub_wls_odometry;

ros::Publisher pub_gnss_raw;
ros::Publisher pub_velocity_from_doppler;

GNSS_Tools m_GNSS_Tools; // utilities

extern void pntposRegisterPub(ros::NodeHandle &n)
{
    // pub_pntpos_odometry = n.advertise<nav_msgs::Odometry>("WLSENURTKLIB", 1000);
    // pub_gnss_raw = n.advertise<nlosExclusion::GNSS_Raw_Array>("GNSSPsrCarRov1", 1000);
    // pub_wls_odometry = n.advertise<nav_msgs::Odometry>("WLSENUGoGPS", 1000);
    pub_velocity_from_doppler = n.advertise<nav_msgs::Odometry>("GNSSDopVelRov1", 1000); // velocity_from_doppler
}

// add by Yixin
/* doppler residuals for LEO obtained from csv ---------------------------------------------------
* compute doppler residuals for LEO satellites
* args   : double *doppler_shifts I   doppler shifts for each satellite (Hz)
*          double *lambdas        I   wavelengths for each satellite (m)
*          int    n                            I   number of observation data
*          double *rs                          I   satellite positions and velocities (ECEF) (m|m/s)
*          double *dts                         I   satellite clock biases and drifts (s)
*          double *rr                          I   receiver position (ECEF) (m)
*          double *x                           I   state vector (receiver velocity and clock bias)
*          double *azel                        I   azimuth/elevation angles (rad)
*          int    *vsat                        I   valid satellite flags
*          double *v                           O   doppler residuals
*          double *H                           O   design matrix
* return : int                                 number of valid observations used in the computation
* notes  : this function does not use obs and nav structures, but directly uses doppler shifts and lambdas
*/
static int resdop_LEO(const double *doppler_shifts, const double *lambdas,
                      int n, const double *rs, const double *dts, const double *rr, const double *x,
                      const double *azel, const int *vsat, double *v, double *H)
{
    double rate, pos[3], E[9], a[3], e[3], vs[3], cosel;
    int i, j, nv = 0;

    trace(3, "resdop_LEO  : n=%d\n", n);

    ecef2pos(rr, pos);
    xyz2enu(pos, E);

    for (i = 0; i < n && i < MAXOBS; i++) {
        double lam = lambdas[i];

        if (doppler_shifts[i] == 0.0 || lam == 0.0 || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0) {
            continue;
        }
        /* line-of-sight vector in ecef */
        cosel = cos(azel[1 + i * 2]);
        a[0] = sin(azel[i * 2]) * cosel;
        a[1] = cos(azel[i * 2]) * cosel;
        a[2] = sin(azel[1 + i * 2]);
        matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);

        /* satellite velocity relative to receiver in ecef */
        for (j = 0; j < 3; j++) vs[j] = rs[j + 3 + i * 6] - x[j];

        /* range rate with earth rotation correction */
        rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * x[0] -
                                                rs[3 + i * 6] * rr[1] - rs[i * 6] * x[1]);

        /* doppler residual */
        v[nv] = -lam * doppler_shifts[i] - (rate + x[3] - CLIGHT * dts[1 + i * 2]);

        /* design matrix */
        for (j = 0; j < 4; j++) H[j + nv * 4] = j < 3 ? -e[j] : 1.0;

        nv++;
    }
    return nv;
}

// add by Yixin
/**
 * @brief Estimate receiver velocity using Doppler shifts for LEO satellites.
 *
 * This function estimates the receiver's velocity based on the provided Doppler shifts,
 * satellite positions, velocities, and other parameters.
 *
 * @param doppler_shifts A vector of Doppler shifts for each satellite, in Hz.
 * @param lambdas A vector of wavelengths for each satellite, in meters.
 * @param n The number of observation data points.
 * @param rs An array of satellite positions and velocities in ECEF coordinates, in meters and meters per second.
 *           The array length should be 6 * n, with each satellite's data in the format [x, y, z, vx, vy, vz].
 * @param dts An array of satellite clock biases and drifts, in seconds.
 *            The array length should be 2 * n, with each satellite's data in the format [dt, dtr].
 * @param rr The receiver's position in ECEF coordinates, in meters. The array length should be 3.
 * @param azel An array of azimuth and elevation angles for each satellite, in radians.
 *             The array length should be 2 * n, with each satellite's data in the format [azimuth, elevation].
 * @param vsat An array of flags indicating whether each satellite is valid. The array length should be n.
 * @param dop_res An array to store the computed Doppler residuals. The array length should be n.
 * @param velocity An array to store the estimated receiver velocity. The array length should be 3.
 */
extern int estvel_LEO(const double *doppler_shifts, const double *lambdas,
                       int n, const double *rs, const double *dts, const double *rr,
                       const double *azel, const int *vsat, double *dop_res, double *velocity)
{
    double x[4] = {0}, dx[4], Q[16], *v, *H;
    int i, j, nv;

    trace(3, "estvel_LEO  : n=%d\n", n);

    v = mat(n, 1);
    H = mat(4, n);

    for (i = 0; i < MAXITR; i++) {

        /* doppler residuals */
        if ((nv = resdop_LEO(doppler_shifts, lambdas, n, rs, dts, rr, x, azel, vsat, v, H)) < 4) {
            break;
        }
        /* least square estimation */
        if (lsq(H, v, 4, nv, dx, Q)) break;

        for (j = 0; j < 4; j++) x[j] += dx[j];

        if (norm(dx, 4) < 1E-6) {
            for (i = 0; i < 3; i++) velocity[i] = x[i];
            break;
        }
    }
    dop_res = v;
    free(v);
    free(H);
}


/**
 * @brief Compute receiver velocity using LEO satellites doppler shift and publish the result.
 *
 * This function computes the receiver's position and velocity based on the provided Doppler shifts,
 * satellite positions, velocities, and other parameters, and publishes the result as a nav_msgs::Odometry message.
 *
 * @param doppler_shifts An array of Doppler shifts for each satellite, in Hz.
 * @param lambdas An array of wavelengths for each satellite, in meters.
 * @param n The number of observation data points.
 * @param rs An array of satellite positions and velocities in ECEF coordinates, in meters and meters per second.
 *           The array length should be 6 * n, with each satellite's data in the format [x, y, z, vx, vy, vz].
 * @param dts An array of satellite clock biases and drifts, in seconds.
 *            The array length should be 2 * n, with each satellite's data in the format [dt, dtr].
 * @param rr The receiver's position in ECEF coordinates, in meters. The array length should be 3.
 * @param azel An array of azimuth and elevation angles for each satellite, in radians.
 *             The array length should be 2 * n, with each satellite's data in the format [azimuth, elevation].
 * @param vsat An array of flags indicating whether each satellite is valid. The array length should be n.
 * @return int Status of the computation (1: success, 0: failure).
 */
extern int pntpos_LEO(const double* doppler_shifts, const double* lambdas,
                      int n, const double *rs, const double *dts, const double *rr,
                      const double *azel, const int *vsat)
{
    double dop_res[n];
    double velocity[3] = {0};

    // Estimate receiver velocity using Doppler shifts for LEO satellites
    estvel_LEO(doppler_shifts, lambdas, n, rs, dts, rr, azel, vsat, dop_res, velocity);

    // Convert ECEF position to ENU position
    // Eigen::Matrix<double, 3, 1> ENU_ref;
    // ENU_ref << ref_lon, ref_lat, ref_alt;
    // Eigen::Matrix<double, 3, 1> ENU;
    // Eigen::Matrix<double, 3, 1> ECEF;
    // ECEF << rr[0], rr[1], rr[2];
    // ENU = m_GNSS_Tools.ecef2enu(ENU_ref, ECEF);

    // Publish the result as nav_msgs::Odometry
    nav_msgs::Odometry odometry;
    odometry.header.frame_id = "map";
    odometry.child_frame_id = "map";
    odometry.pose.pose.position.x = 0;//ENU(0);
    odometry.pose.pose.position.y = 0;//ENU(1);
    odometry.pose.pose.position.z = 0;//ENU(2);
    odometry.twist.twist.linear.x = velocity[0];
    odometry.twist.twist.linear.y = velocity[1];
    odometry.twist.twist.linear.z = velocity[2];

    // Optionally, set the covariance of the twist (velocity) based on the Doppler residuals
    odometry.twist.covariance[0] = norm(dop_res, n);

    pub_pntpos_odometry.publish(odometry);

    return 1; // Return success status
}
