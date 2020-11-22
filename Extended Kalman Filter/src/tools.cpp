#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    // Check the validity of the following inputs:
    // The estimation vector size should not be zero
    if(estimations.size() == 0){
        cout << "Input is empty" << endl;
        return rmse;
    }
    // The estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()){
        cout << "Invalid estimation or ground_truth. Data should have the same size" << endl;
        return rmse;
    }
    // Accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        VectorXd residual = (estimations[i] - ground_truth[i]);
        // Coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    // Calculate the mean
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
    
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
 
    //Initiallize the Jacobian Matrix
    MatrixXd Jacobian = MatrixXd::Zero(3,4);
    
    //Initialize state parameters px, py, ux, uy
    float px = x_state(0);
    float py = x_state(1);
    float ux = x_state(2);
    float uy = x_state(3);
    
    
    //calculate repeated variables
    float squared_sum = px*px + py*py;
    float root_of_squared_sum = sqrt(squared_sum);
    float three_out_of_two_sum = squared_sum * root_of_squared_sum;
    
    //check if it is close to zero
    if (squared_sum>0.000001){
        Jacobian<<px/root_of_squared_sum, py/root_of_squared_sum, 0, 0,
                -py/squared_sum, px/squared_sum, 0, 0,
        py*(ux*py-uy*px)/three_out_of_two_sum, px*(uy*px-ux*py)/three_out_of_two_sum, px/root_of_squared_sum, py/root_of_squared_sum;
    }
    
    return Jacobian;
    
}
