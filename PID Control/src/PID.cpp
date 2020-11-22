#include "PID.h"
#include <math.h>
#include <vector>
#include <limits>
#include <iostream>

using namespace std;

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd) {
    
    this->Kp=Kp;
    this->Ki=Ki;
    this->Kd=Kd;
    p_error = 0.0;
    i_error = 0.0;
    d_error = 0.0;
    //best_error = 100000000; //Really large value -- bad practise but c++
    best_error = std::numeric_limits<int>::max();
    sum_cte = 0;
    prev_cte = 0;
    
}

void PID::UpdateError(double cte) {
    
    sum_cte += cte;
    p_error = -1 * Kp * cte;
    i_error = -1 * Ki * sum_cte;
    d_error = -1 * Kd * (cte - prev_cte);
    prev_cte = cte;

}

double PID::TotalError() {
    //Implement the PID formula
    return p_error + i_error + d_error;
}


//Defined a new function Twiddle as taught in the classroom - p and dp were chosen to have the same form as the algorithm in the classroom.
// p={Kp, Kd, Ki} and dp={their errors}
void PID::Twiddle(double cte)
{
    //Define the p and dp vectors and tolerance
    std::vector<double> p = {this->Kp,this->Kd,this->Ki};
    std::vector<double> dp = {1, 1, 1};
    double tolerance =  5;
    
    //Calculate sum
    double sum = 0 ;
    for (int i = 0; i < dp.size(); i++){
        sum += dp[i];
    }
    
    //Implement Twiddle algorithm
   while (sum > tolerance){
        for (int i=0; i<p.size(); i++){
            p[i] = p[i] + dp[i];
            if (fabs(cte)<fabs(best_error)){
                best_error = cte;
                dp[i] = dp[i] * 1.1;
            }else{
                p[i] = p[i] - 2 * dp[i];
                if (fabs(cte)<fabs(best_error)){
                    best_error = cte;
                    dp[i] = dp[i] * 1.1;
                }else{
                    p[i] = p[i] + dp[i];
                    dp[i] = dp[i] * 0.9;
                }
            }
        }
       
       //Calculate sum
        sum = 0 ;
        for (int i = 0; i < dp.size(); i++){
            sum += dp[i];
        }
    }
   
    //Update the parameters
    Kp = p[0];
    Kd = p[1];
    Ki = p[2];
    p_error = dp[0];
    d_error = dp[1];
    i_error = dp[2];
}


