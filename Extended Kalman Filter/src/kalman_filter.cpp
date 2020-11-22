#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    
    //Predict with the equaltions x'=Fx and P'=FPF(t)+Q
    x_ = F_ * x_ ;
    MatrixXd F_T = F_.transpose();
    P_ =  F_ * P_ * F_T + Q_ ;
}

void KalmanFilter::Update(const VectorXd &z) {
  
    //Update using the following equations:
    /*
     y=z-Hx'
     S=HP'H(t)+R
     K=P'H(t)S^(-1)
     x=x'+Ky
     P=(I-KH)P'
     */
    VectorXd y = z - H_ * x_;
    MatrixXd H_T = H_.transpose();
    MatrixXd S = H_ * P_ * H_T + R_ ;
    MatrixXd K = P_ * H_T * S.inverse();
    x_ = x_ + K * y ;
    int size = x_.size();
    MatrixXd I = MatrixXd::Identity (size,size);
    P_ = (I - K * H_) * P_ ;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
    //Update using the same equations as above. The difference here is that we need to calculated the polar coordinates in order to find the y vector since z contains the polar coordinates.
    
    //Variables initiallization
    float px = x_[0];
    float py = x_[1];
    float ux = x_[2];
    float uy = x_[3];
    float rho = sqrt(pow(px,2)+pow(py,2));
    float phi = atan2(py,px);
    float rhodot = (ux * px + uy * py) / rho;
    
    VectorXd z_prediction (3);
    z_prediction << rho, phi, rhodot;
    VectorXd y = z - z_prediction;
    
    //Normalization is requires here When calculating phi in y = z - h(x) for radar measurements, the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi. The Kalman filter is expecting small angle values between the range -pi and pi. When working in radians, we can add 2π or subtract 2π until the angle is within the desired range.
    
    float double_pi = M_PI*2;
    
    while (y[1] > M_PI){
        y[1] = y[1] - double_pi;
    }
    
    while (y[1] < - M_PI){
        y[1] = y[1] + double_pi;
    }
    
    //Same equation as above
    MatrixXd H_T = H_.transpose();
    MatrixXd S = H_ * P_ * H_T + R_ ;
    MatrixXd K = P_ * H_T * S.inverse();
    x_ = x_ + K * y ;
    int size = x_.size();
    MatrixXd I = MatrixXd::Identity (size,size);
    P_ = (I - K * H_) * P_ ;
    
    
}
