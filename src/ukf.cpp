#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
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
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   n_x_ = 5;
   n_aug_ = 7;
   n_z_radar = 3;
   n_z_lidar = 2;
   lambda_ = 3 - n_aug_;
   
  // create augmented mean vector 
  x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  // create and initialise Weight vector
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  weights_(0) = lambda_ / (lambda_ + n_aug_ );
  for (int i=1 ; i< (2 * n_aug_ + 1) ; i++)
    weights_(i) = 0.5 / (lambda_ + n_aug_);
   
   
}

UKF::~UKF() {}

void UKF::GenerateSigmaPoint(){
  // create augmented mean state
  
  x_aug.head(n_x_) = x_;
  // add noise Va / Vaa
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;
  
  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) =  std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
  
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  
  Xsig_aug.col(0) = x_aug;
  //std::cout << "x_aug = " << std::endl << Xsig_aug << std::endl;
  
  for(int i=0;i < n_aug_;i++){
      
      Xsig_aug.col(i+1)           = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  
	
}


void UKF::AugmentSigmaPoint(){
}

void UKF::PredictSigmaPoint(){
}

void UKF::PredictMeanCovariance(){
}





void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
	  float rho_mea=meas_pack.raw_measurements_[0];
	  float theta_mea=meas_pack.raw_measurements_[1];
	  float rhodot_mea=meas_pack.raw_measurements_[2];
	  
	  float px = rho_mea * cos( theta_mea);
	  float py = rho_mea * sin(theta_mea);
	  float vx = rhodot_mea * cos( theta_mea);
	  float vy = rhodot_mea * cos( theta_mea);
	  
	  
	  x_ << px, 
		    py, 
            sqrt(vx * vx + vy * vy), 
            0,
			0; 

    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
	  x_ << meas_pack.raw_measurements_[0], 
            meas_pack.raw_measurements_[1], 
            0, 
            0,
			0;

    }
	
	
	if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
	VectorXd z=VectorXd(3);
	z<< meas_pack.raw_measurements_[0],meas_pack.raw_measurements_[1],meas_pack.raw_measurements_[2];
	

  } else {
    // TODO: Laser updates
	VectorXd z=VectorXd(2);
	z<< meas_pack.raw_measurements_[0],meas_pack.raw_measurements_[1];
	

  }
   
   
   
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
   // Predict Mean
   for (int i=0 ; i< n_x_ ; i++)
    x_(i) = Xsig_pred.row(i) * weights_;
  

   // predict state covariance matrix
   for (int i=0 ; i< n_x_ ; i++){
     VectorXd term1 = Xsig_pred.row(i) - x_;
     P_.row(i) = term1 * term1.transpose() * weights_;
   }
   
   
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
   for (int i=0; i<2 * n_aug_ + 1; ++i){
      
      Zsig(0,i) = Xsig_pred(0,i);
      Zsig(1,i) = Xsig_pred(1,i);

  }
  
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   
  // transform sigma points into measurement space
  double rho,phi,rho_dot;
  double x_px, x_py, x_vel, x_phi, x_phi_dot;
  
  for (int i=0; i<2 * n_aug_ + 1; ++i){
      
      x_px = Xsig_pred(0,i);
      x_py = Xsig_pred(1,i);
      x_vel = Xsig_pred(2,i);
      x_phi = Xsig_pred(3,i);
      x_phi_dot = Xsig_pred(4,i);
      
      
      rho = sqrt(x_px * x_px + x_py * x_py );
      phi = atan(x_py/x_px);
      rho_dot = (x_px * cos(x_phi) * x_vel + x_py * sin(x_phi) * x_vel) / rho;
      
      Zsig(0,i) = rho;
      Zsig(1,i) = phi;
      Zsig(2,i) = rho_dot;
      
  }
   
  // calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar);
  z_pred.fill(0.0);
  for (int j=0; j< 2 * n_aug_ + 1; ++j){
      z_pred = z_pred + Zsig.col(j) * weights_(j);
      
  }
  
  // calculate innovation covariance matrix S
  
  MatrixXd R = MatrixXd(n_z_radar,n_z_radar);
  R<< std_radr_ * std_radr_,0,0,
      0,std_radphi_ * std_radphi_,0,
      0,0,std_radrd_ * std_radrd_;
  
  MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
  S.fill(0.0);      
  for (int k=0; k < 2 * n_aug_ + 1; ++k){

      
      VectorXd diff = Zsig.col(k) - z_pred;
      S = S + diff * diff.transpose() * weights_(k);

  }
  S = S + R;  
   
   
   
   
   
}