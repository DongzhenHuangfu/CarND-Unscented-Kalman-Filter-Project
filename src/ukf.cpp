#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // State dimension
  int n_x_ = 5;

  // Augmented state dimension
  int n_aug_ = 7;

  // Sigma point spreading parameter
  double lambda_ = 3 - a_aug_;

  //create vector for weights
  VectorXd weights_ = VectorXd(2*n_aug+1);

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_){
    cout << "UKF: " << endl;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float theta_d = meas_package.raw_measurements_(2);

      float px = rho * cos(theta);
      float py = rho * sin(theta);

      x_ << px, py, 0, 0, 0;
    }

    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){

      float px = meas_package.raw_measurements_(0);
      float py = meas_package.raw_measurements_(1);

      x_ << px, py, 0, 0, 0;
    }
    time_us_ = meas_package.timestamp_;

    P << 0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
        -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
         0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
        -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
        -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //set weights
  double weights0 = lambda/(lambda + n_aug);
  weights_(0) = weights0;
  x_.fill(0.0);
  P_.fill(0.0);
  for(int i = 1; i < 2*n_aug+1; i++)
  {
      double weight = 0.5/(lambda + n_aug);
      weights_(i) = weight;
  }

    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  cout<<"Prediction:"<<endl;
  cout<<"x_: "<<x_<<endl;

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateLidar(meas_package);
    return;
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
    return;
  }
}

/* function for generating the sigmoid points*/

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){

  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a * std_a;
  P_aug(6,6) = std_yawdd * std_yawdd;

  MatrixXd A_aug = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug; i++)
  {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug) * A_aug.col(i);
      Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * A_aug.col(i);
  }

  *Xsig_out = Xsig_aug;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // generate the sigmoid points

  MatrixXd x_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(*x_aug);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  for(int i = 0; i < 2*n_aug+1; i++)
  {
      float px = Xsig_aug(0,i);
      float py = Xsig_aug(1,i);
      float v = Xsig_aug(2,i);
      float yaw = Xsig_aug(3,i);
      float yaw_d = Xsig_aug(4,i);
      float mu_a = Xsig_aug(5, i);
      float mu_yawdd = Xsig_aug(6,i);
      
      VectorXd d1(5);
      VectorXd d2(5); 
      if(fabs(yaw_d) < 0.0001)
      {
          d1 << v * cos(yaw) * delta_t,
                v * sin(yaw) * delta_t,
                0.0,
                0.0,
                0.0;
      }
      else
      {
          d1 << v/yaw_d*(sin(yaw+yaw_d*delta_t)-sin(yaw)),
                v/yaw_d*(-cos(yaw+yaw_d*delta_t) + cos(yaw)),
                0.0,
                yaw_d * delta_t,
                0.0;
      }
      d2 << 0.5*delta_t*delta_t*cos(yaw)*mu_a,
            0.5*delta_t*delta_t*sin(yaw)*mu_a,
            delta_t * mu_a,
            0.5 * delta_t*delta_t*mu_yawdd,
            delta_t*mu_yawdd;
        
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + d1 + d2; 
      
  }

  //predict state mean
  for(int i = 0; i < 2*n_aug+1;i++)
  {
      x_ += weights_(i)*Xsig_pred.col(i);
  }
  //predict state covariance matrix
  for(int i = 0; i < 2*n_aug+1; i++)
  {
      VectorXd x_diff = Xsig_pred.col(i)-x;
      while(x_diff(3)<= -atan(1)*4)
      {
          x_diff(3) += 8*atan(1);
      }
      while(x_diff(3) > atan(1)*4)
      {
          x_diff(3) -= 8*atan(1);
      }
      P_ += weights_(i)*x_diff*x_diff.transpose();
  }

}


void UKF::UpdateUKF(int n_z, MatrixXd Zsig, MatrixXd z_pred, VectorXd z){

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(x_diff(3) <= -4*atan(1)) x_diff(3) += 8*atan(1);
      while(x_diff(3) > 4*atan(1)) x_diff(3) -= 8*atan(1);
      while(z_diff(1) <= -4*atan(1)) z_diff(1) += 8*atan(1);
      while(z_diff(1) > 4*atan(1)) z_diff(1) -= 8*atan(1);
      Tc += weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z);

  MatrixXd S = MatrixXd(n_z,n_z);

  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      double px = x_(0,i);
      double py = x_(1,i);
      
      if(fabs(rho)<0.001)
      {
          rho_d = (v*px*cos(yaw) + v*py*sin(yaw))/0.001;
      }
      else
      {
          rho_d = (v*px*cos(yaw) + v*py*sin(yaw))/rho;
      }
      Zsig.col(i) << px,
                     py;
                     
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(z_diff(1)<=-4*atan(1))
      {
          z_diff(1) = z_diff(1) + 8*atan(1);
      }
      while(z_diff(1)>4*atan(1))
      {
          z_diff(1) = z_diff(1) - 8*atan(1);
      }
      S = S + weights(i)*z_diff*z_diff.transpose();
  }
  
  MatrixXd R(n_z, n_z);
  R << std_laspx_*std_laspx_, 0.0,
       0.0, std_laspy_*std_laspy_,

  S = S + R;

  VectorXd z = VerctorXd(3);

  double px = meas_package.raw_measurements(0);
  double py = meas_package.raw_measurements(1);
  z << px,
       py;

  float e = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  NIS_LIDAR.push_back(e);

  UpdateUKF(n_z, Zsig, z_pred, z);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z);

  MatrixXd S = MatrixXd(n_z,n_z);

  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      double px = x_(0,i);
      double py = x_(1,i);
      double v = x_(2,i);
      double yaw = x_(3,i);
      double yaw_d = x_(4,i);
      
      double rho = sqrt(px*px+py*py);
      double theta = atan2(py,px);
      double rho_d;
      
      if(fabs(rho)<0.001)
      {
          rho_d = (v*px*cos(yaw) + v*py*sin(yaw))/0.001;
      }
      else
      {
          rho_d = (v*px*cos(yaw) + v*py*sin(yaw))/rho;
      }
      Zsig.col(i) << rho,
                     theta,
                     rho_d;
                     
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(z_diff(1)<=-4*atan(1))
      {
          z_diff(1) = z_diff(1) + 8*atan(1);
      }
      while(z_diff(1)>4*atan(1))
      {
          z_diff(1) = z_diff(1) - 8*atan(1);
      }
      S = S + weights(i)*z_diff*z_diff.transpose();
  }
  
  MatrixXd R(n_z, n_z);
  R << std_radr_*std_radr_, 0.0, 0.0,
       0.0, std_radphi_*std_radphi_, 0.0,
       0.0, 0.0, std_radrd_*std_radrd_;
  S = S + R;

  VectorXd z = VerctorXd(3);

  double rho = meas_package.raw_measurements(1);
  double phi = meas_package.raw_measurements(2);
  double rho_dot = meas_package.raw_measurements(3);
  z << rho,
       phi,
       rho_dot;

  float e = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  NIS_RADAR.push_back(e);

  UpdateUKF(n_z, Zsig, z_pred, z);

  return;
}
