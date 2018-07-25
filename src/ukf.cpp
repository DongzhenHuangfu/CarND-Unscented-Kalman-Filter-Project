#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream> 

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

  // State dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  count_Lidar_ = 0;

  count_Lidar_up_ = 0;

  count_Radar_ = 0;

  count_Radar_up_ = 0;
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
  // cout<<"in ukf process!"<<endl;
  if(!is_initialized_){
    // cout << "UKF: " << endl;

    weights_.fill(0.0);
    x_.fill(0.0);

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      // float theta_d = meas_package.raw_measurements_(2);

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

  //set weights
  double weights0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weights0;
  for(int i = 1; i < 2*n_aug_+1; i++)
  {
      double weight = 0.5/(lambda_ + n_aug_);
      weights_(i) = weight;
  }

  is_initialized_ = true;
    
  return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // cout<<"before prediction!"<<endl;

  Prediction(dt);

  // cout<<"Prediction:"<<endl;
  // cout<<"x_: "<<endl<<x_<<endl<<endl;

  if((meas_package.sensor_type_ == MeasurementPackage::RADAR) & (use_radar_)){
    UpdateRadar(meas_package);
    return;
  }
  else if((meas_package.sensor_type_ == MeasurementPackage::LASER) & (use_laser_)){
    UpdateLidar(meas_package);
    return;
  }
}

/* function for generating the sigmoid points*/

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){

  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd A_aug = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * A_aug.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A_aug.col(i);
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
  // cout<<"0"<<endl;

  MatrixXd x_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  AugmentedSigmaPoints(&x_aug);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      float px = x_aug(0,i);
      float py = x_aug(1,i);
      float v = x_aug(2,i);
      float yaw = x_aug(3,i);
      float yaw_d = x_aug(4,i);
      float mu_a = x_aug(5, i);
      float mu_yawdd = x_aug(6,i);

      // predict state values
      double px_p, py_p;

      if(fabs(yaw_d) < 0.0001)
      {
        px_p = px + v * delta_t * cos(yaw);
        py_p = py + v * delta_t * sin(yaw);
      }
      else
      {
        px_p = px + v / yaw_d * ( sin(yaw + yaw_d * delta_t) - sin(yaw));
        py_p = py + v / yaw_d * ( cos(yaw) - cos(yaw + yaw_d * delta_t));
      }

      double v_p = v;
      double yaw_p = yaw + yaw_d * delta_t;
      double yawd_p = yaw_d;
        
      // add noise
      px_p += 0.5 * mu_a * delta_t * delta_t * cos(yaw);
      py_p += 0.5 * mu_a * delta_t * delta_t * sin(yaw);
      v_p += mu_a *delta_t;

      yaw_p += 0.5 * mu_yawdd * delta_t * delta_t;
      yawd_p += mu_yawdd * delta_t;

      // write the prediction
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
  }

  // cout<<"Xsig_pred_:"<<endl<<Xsig_pred_<<endl;

  // cout<<"1"<<endl;
  x_.fill(0.0);
  P_.fill(0.0);

  //predict state mean
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      x_ += weights_(i)*Xsig_pred_.col(i);
  }

  // cout<<"2"<<endl;
  //predict state covariance matrix
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      // cout<<"x_diff:"<<x_diff<<endl;

      while(x_diff(3) <= -atan(1)*4)
      {
          x_diff(3) += 8*atan(1);
      }
      // cout<<"2in1-"<<i<<endl;
      while(x_diff(3) > atan(1)*4)
      {
          x_diff(3) -= 8*atan(1);
      }
      // cout<<"2in2-"<<i<<endl;
      P_ += weights_(i)*x_diff*x_diff.transpose();
  }
  // cout<<"3"<<endl;
  return;
}


void UKF::UpdateUKF(int n_z, MatrixXd Zsig, MatrixXd z_pred, VectorXd z, MatrixXd S){

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for(int i = 0; i< 2 * n_aug_ + 1; i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(x_diff(3) <= -4*atan(1)) x_diff(3) += 8*atan(1);
      while(x_diff(3) > 4*atan(1)) x_diff(3) -= 8*atan(1);
      while(z_diff(1) <= -4*atan(1)) z_diff(1) += 8*atan(1);
      while(z_diff(1) > 4*atan(1)) z_diff(1) -= 8*atan(1);
      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();
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

  MatrixXd S = MatrixXd(n_z, n_z);

  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);

      Zsig.col(i) << px,
                     py;
                     
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // cout<<"z Prediction of Lidar:"<<endl<<z_pred<<endl<<endl;

  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;

      S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z, n_z);

  R << std_laspx_*std_laspx_, 0.0,
       0.0, std_laspy_*std_laspy_;

  S += R;

  VectorXd z = VectorXd(n_z);

  double px = meas_package.raw_measurements_(0);
  double py = meas_package.raw_measurements_(1);
  z << px,
       py;

  float e = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  count_Lidar_++;
  if(e >= 3){
    count_Lidar_up_++;
  }
  float percent = (float)count_Lidar_up_ / count_Lidar_;

  cout<<"Lidar NIS, persent of over 90%: "<<percent<<", e:"<<e<<endl;
  cout<<"number of total:"<<count_Lidar_<<", number of up:"<<count_Lidar_up_<<endl;

  UpdateUKF(n_z, Zsig, z_pred, z, S);

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
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
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

  // cout<<"z Prediction of Radar:"<<endl<<z_pred<<endl<<endl;

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
      S = S + weights_(i)*z_diff*z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0.0, 0.0,
       0.0, std_radphi_*std_radphi_, 0.0,
       0.0, 0.0, std_radrd_*std_radrd_;
  S += R;

  VectorXd z = VectorXd(n_z);

  double rho = meas_package.raw_measurements_(0);
  double phi = meas_package.raw_measurements_(1);
  double rho_dot = meas_package.raw_measurements_(2);
  z << rho,
       phi,
       rho_dot;

  float e = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  count_Radar_++;
  if(e >= 7.8){
    count_Radar_up_++;
  }
  float percent = (float)count_Radar_up_ / count_Radar_;

  cout<<"Radar NIS, persent of over 90%: "<<percent<<", e:"<<e<<endl;
  cout<<"number of total:"<<count_Radar_<<", number of up:"<<count_Radar_up_<<endl;

  UpdateUKF(n_z, Zsig, z_pred, z, S);

  return;
}
