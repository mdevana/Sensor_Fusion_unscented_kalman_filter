#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	
	
  // if this is false, then state vector is initialised using measurements	
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  
  P_ = MatrixXd::Identity(5,5);

  // create sigma point matrix
  Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  // create and initialise Weight vector
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  weights_(0) = lambda_ / (lambda_ + n_aug_ );
  weights_.tail(2*n_aug_).fill(0.5 / (n_aug_ + lambda_));
  /*for (int i=1 ; i< (2 * n_aug_ + 1) ; i++)
    weights_(i) = 0.5 / (lambda_ + n_aug_);*/

  H_laser_ = MatrixXd(n_z_lidar, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
             0, 1, 0, 0, 0;
			 
  R_laser_ = MatrixXd(n_z_lidar,n_z_lidar);
  R_laser_<< std_laspx_ * std_laspx_,0,
      0,std_laspy_ * std_laspy_;
	  
  R_radar_ = MatrixXd(n_z_radar,n_z_radar);
  R_radar_<< std_radr_ * std_radr_,0,0,
      0,std_radphi_ * std_radphi_,0,
      0,0,std_radrd_ * std_radrd_;
	  
  
   
   
}

UKF::~UKF() {}

void UKF::Prediction(double dt){
	
  
   AugmentSigmaPoint();
   PredictSigmaPoint(dt);
   PredictMeanCovariance();
	
	
}

void UKF::AugmentSigmaPoint(){
  // create augmented mean state
  
  x_aug.head(n_x_) = x_;
  // add noise Va / Vaa
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;
  
  
  
  // create augmented covariance matrix
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) =  std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
  
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  
  Xsig_aug.col(0) = x_aug;
  
  
  for(int i=0;i < n_aug_;i++){
      
      Xsig_aug.col(i+1)            = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  
	
}


void UKF::PredictSigmaPoint(double delta_t){

  Xsig_pred_.fill(0);
  double px,py,v,phi,phidot,std_a, std_yaw;
  double px_i, py_i,v_i, phi_i, phidot_i;
  
  
  for(int i=0;i < 2 * n_aug_ + 1;i++){
      px = Xsig_aug(0,i);
      py = Xsig_aug(1,i);
      v = Xsig_aug(2,i);
      phi = Xsig_aug(3,i);
      phidot = Xsig_aug(4,i);
      std_a = Xsig_aug(5,i);
      std_yaw = Xsig_aug(6,i);
      
      if (fabs(phidot) > 0.001){
          px_i = v/phidot * (sin( phi + phidot * delta_t) - sin(phi));
          py_i = v/phidot * (-1*cos( phi + phidot * delta_t) + cos(phi));
          v_i = 0;
          phi_i = phidot *delta_t;
          phidot_i = 0;
          
          
      }
      else{
           px_i = v * cos(phi) * delta_t;
           py_i = v * sin(phi) * delta_t;
           v_i = 0;
           phi_i = phidot * delta_t;
           phidot_i = 0;
           
      }
      
      px_i += px + (0.5 * (delta_t * delta_t) * cos(phi) * std_a); 
      py_i += py + (0.5 * (delta_t * delta_t) * sin(phi) * std_a); 
      v_i  += v  + (1.0 * (delta_t) * std_a);
      phi_i += phi + (0.5 * (delta_t * delta_t) * std_yaw);
      phidot_i += phidot + (1.0 * (delta_t) * std_yaw);
      
      Xsig_pred_(0,i) = px_i;
      Xsig_pred_(1,i) = py_i;
      Xsig_pred_(2,i) = v_i;
      Xsig_pred_(3,i) = phi_i;
      Xsig_pred_(4,i) = phidot_i;

  }
	
	
	
}



void UKF::PredictMeanCovariance(){
	/**
   * Estimate the object's location. 
   * Modify the state vector, x_ and
   * and the state covariance matrix.
   */
  
  // Predict mean of state vector based on Xsig_pred
  
  x_.fill(0);
  for (int i=0 ; i< 2 * n_aug_ + 1 ; i++)
    x_ = x_ + Xsig_pred_.col(i) * weights_(i);
  

  // predict state covariance matrix
  P_.fill(0);
  for (int i=0 ; i< 2 * n_aug_ + 1; i++){
    VectorXd term1 = Xsig_pred_.col(i) - x_;
    term1(3) = WrapAngle(term1(3));
    P_ = P_ + weights_(i) * term1 * term1.transpose();

  }

  
}

double UKF::WrapAngle(double angleValue){
	
	while (angleValue > M_PI ) angleValue-= 2.0 * M_PI;
	while (angleValue < -1 * M_PI ) angleValue+= 2.0 * M_PI;
	
	return (angleValue);
	
}




void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
   * switch between lidar and radar measurements.
   */
   

   
   if (!is_initialized_) {
	// Initialise the state vector   
	if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      // and initialize state.
	  
	  float rho_mea=meas_pack.raw_measurements_[0];
	  float theta_mea=meas_pack.raw_measurements_[1];
	  float rhodot_mea=meas_pack.raw_measurements_[2];
	  
	  
	  
	  float px = rho_mea * cos(theta_mea);
	  float py = rho_mea * sin(theta_mea);
	  float vx = rhodot_mea * cos( theta_mea);
	  float vy = rhodot_mea * cos( theta_mea);
	  float v = sqrt(vx * vx + vy * vy);
	  	  
	  x_ << px, 
		    py, 
            v, 
            0,
			0;
	  P_<< std_radr_ * std_radr_,0,0,0,0,
           0,std_radr_ * std_radr_,0,0,0,
           0,0,std_radrd_ * std_radrd_,0,0,
		   0,0,0,std_radphi_ * std_radphi_,0,
		   0,0,0,0,std_radphi_ * std_radphi_;
	is_initialized_ = true;

    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
	  x_ << meas_pack.raw_measurements_[0], 
            meas_pack.raw_measurements_[1], 
            0, 
            0,
			0;

	  P_<< std_laspx_ * std_laspx_,0,0,0,0,
           0, std_laspy_  *  std_laspy_ ,0,0,0,
           0,0,1,0,0,
		   0,0,0,1,0,
		   0,0,0,0,1;
    is_initialized_ = true; 
    }
	
	previous_timestamp_=meas_pack.timestamp_;
	
	
   }
   
   
   // Stage after initialisation
   else {

   float dt = (meas_pack.timestamp_-previous_timestamp_)/1000000.0;
   previous_timestamp_=meas_pack.timestamp_;
   Prediction(dt);
   
	
	if ((meas_pack.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_ == true)) {
    // Radar updates
		UpdateRadar(meas_pack);
   } else if ((use_laser_) && (meas_pack.sensor_type_ == MeasurementPackage::LASER)){

		//UpdateLidar(meas_pack);
		UpdateLidar_linear(meas_pack); // Use linear transformation algortihm for Laser

   }

 }// End Else module
   
   
   
}
void UKF::UpdateLidar_linear(MeasurementPackage meas_pack){
  
	  VectorXd z=VectorXd(n_z_lidar);
      z<<meas_pack.raw_measurements_[0],meas_pack.raw_measurements_[1];
				  
	  VectorXd y = z - H_laser_ * x_ ;
	  
      MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_ ;
	  
      MatrixXd K = P_ * H_laser_.transpose() * S.inverse();
   
      long x_size = x_.size();
      MatrixXd I = MatrixXd::Identity(x_size, x_size);
      
      x_ = x_ + ( K * y) ;
      P_ = (I - K * H_laser_) * P_;
				  
}



void UKF::UpdateLidar(MeasurementPackage meas_pack) {
  /**
   *  Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
   
   
   MatrixXd Zsig = MatrixXd(n_z_lidar, 2 * n_aug_ + 1);
   for (int i=0; i<2 * n_aug_ + 1; ++i){
      
      Zsig(0,i) = Xsig_pred_(0,i);
      Zsig(1,i) = Xsig_pred_(1,i);

  }
  
  // calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_lidar);
  z_pred.fill(0.0);
  for (int j=0; j< 2 * n_aug_ + 1; ++j){
      z_pred = z_pred + Zsig.col(j) * weights_(j);
      
  }
  
  // calculate innovation covariance matrix S
  
  

  
  MatrixXd S = MatrixXd(n_z_lidar,n_z_lidar);
  S.fill(0.0);      
  for (int k=0; k < 2 * n_aug_ + 1; ++k){

      
      VectorXd diff = Zsig.col(k) - z_pred;
      S = S + diff * diff.transpose() * weights_(k);

  }
  S = S + R_laser_;  
  
  // calculate Cross-Correlation matrix TC
  
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar);
  Tc.fill(0);
  
  VectorXd diff_X;
  VectorXd diff_Z;
  
  for(int k = 0; k < 2 * n_aug_ + 1 ; k++  ) {
      diff_X = Xsig_pred_.col(k) - x_;
      diff_Z = Zsig.col(k) - z_pred;
      
      Tc = Tc + weights_(k) * diff_X * diff_Z.transpose();
      
  }
  
  MatrixXd Kgain;
  VectorXd diff_z;
  
  Kgain = Tc * S.inverse();
  
  VectorXd z=VectorXd(n_z_lidar);
  z<< meas_pack.raw_measurements_[0],meas_pack.raw_measurements_[1];
  diff_z = z - z_pred;
  
  x_ = x_ + Kgain * diff_z;
  P_ = P_ - Kgain * S * Kgain.transpose();
  
  NIS_lidar_ = diff_z.transpose() * S.inverse() * diff_z;

}

void UKF::UpdateRadar(MeasurementPackage meas_pack) {
  /**
   * Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // transform sigma points into measurement space
  double rho,phi,rho_dot;
  double x_px, x_py, x_vel, x_phi, x_phi_dot;
  
  MatrixXd Zsig = MatrixXd(n_z_radar, 2 * n_aug_ + 1);
  
  for (int i=0; i<2 * n_aug_ + 1; ++i){
      
      x_px = Xsig_pred_(0,i);
      x_py = Xsig_pred_(1,i);
      x_vel = Xsig_pred_(2,i);
      x_phi = Xsig_pred_(3,i);
      x_phi_dot = Xsig_pred_(4,i);
      
      
      rho = sqrt(x_px * x_px + x_py * x_py );
      phi = atan2(x_py,x_px);
	  
	  if (rho < 0.001)
			rho_dot = (x_px * cos(x_phi) * x_vel + x_py * sin(x_phi) * x_vel) / 0.001;
      else		  
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
  
  MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
  S.fill(0.0);      
  for (int k=0; k < 2 * n_aug_ + 1; ++k){

      
      VectorXd diff = Zsig.col(k) - z_pred;
	  diff(1)=WrapAngle(diff(1));
	  
  
      S = S + diff * diff.transpose() * weights_(k);

  }
  S = S + R_radar_;  
  
  // calculate Cross-Correlation matrix TC
  
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);
  Tc.fill(0);
  
  VectorXd diff_X;
  VectorXd diff_Z;
  
  for(int k = 0; k < 2 * n_aug_ + 1 ; k++  ) {
      diff_X = Xsig_pred_.col(k) - x_;
      diff_Z = Zsig.col(k) - z_pred;
	  
	  diff_Z(1)=WrapAngle(diff_Z(1));
	  
	  diff_X(3)=WrapAngle(diff_X(3));
      
      Tc = Tc + weights_(k) * diff_X * diff_Z.transpose();
      
  }
  
  MatrixXd Kgain;
  VectorXd diff_z;
  
  Kgain = Tc * S.inverse();
  
  VectorXd z=VectorXd(3);
  z<< meas_pack.raw_measurements_[0],meas_pack.raw_measurements_[1],meas_pack.raw_measurements_[2];
  
  diff_z = z - z_pred;  
  diff_z(1)=WrapAngle(diff_z(1));
  
  x_ = x_ + Kgain * diff_z;
  P_ = P_ - Kgain * S * Kgain.transpose();
  
  NIS_radar_ = diff_z.transpose() * S.inverse() * diff_z;
  

  //std::cout << "Theta value = " <<theta_mea * 180 / M_PI <<" ";
  //std::cout<<x_(0)<<" "<<x_(1)<<" "<< x_(3)<< std::endl;
  //std::cout << "p_ = " << std::endl << P_ << std::endl;

}