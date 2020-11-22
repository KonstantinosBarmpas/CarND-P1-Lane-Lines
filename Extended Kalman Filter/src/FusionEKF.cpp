#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  //measurement noise covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement noise covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    
    //initialize measurement function H - laser
    H_laser_<<1, 0, 0, 0,
        0, 1, 0, 0;

    //initialize measurement function H - radar
    Hj_ = MatrixXd::Zero(3,4);
    
    //initialize the F matrix
    ekf_.F_<< 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
    
    //initialize P matrix
    ekf_.P_ = MatrixXd::Identity(4,4);
    
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        float rho_measured = measurement_pack.raw_measurements_(0);
        float phi_measured = measurement_pack.raw_measurements_(1);
        float rhodot_measured = measurement_pack.raw_measurements_(2); //I didint use it somewhere
        ekf_.x_ << rho_measured*cos(phi_measured), rho_measured*sin(phi_measured), 0, 0; //px, py, ux, uy
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        float px_measured = measurement_pack.raw_measurements_(0);
        float py_measured = measurement_pack.raw_measurements_(1);
        ekf_.x_[0] = px_measured; //px
        ekf_.x_[1] = py_measured; //py
    }

      previous_timestamp_=measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    float dt_micro = (measurement_pack.timestamp_ - previous_timestamp_);
    float dt = dt_micro / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    float noise_ax = 9.0;
    float noise_ay = 9.0;
    
    //Update the transition F matrix
    ekf_.F_<< 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;
    
    //Update the process noise covariance matrix using noise_ax = 9 and noise_ay = 9
    ekf_.Q_<< pow(dt,4.0) / 4.0 * noise_ax, 0, pow(dt,3.0) / 2.0 * noise_ax, 0,
            0, pow(dt,4.0) / 4.0 * noise_ay, 0, pow(dt,3.0) / 2.0 * noise_ay,
            pow(dt,3.0) / 2.0 * noise_ax, 0, pow(dt,2.0) * noise_ax, 0,
            0, pow(dt,3.0) / 2.0 * noise_ay, 0, pow(dt,2.0) * noise_ay;
    

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.H_=tools.CalculateJacobian(ekf_.x_);
      Hj_=ekf_.H_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
