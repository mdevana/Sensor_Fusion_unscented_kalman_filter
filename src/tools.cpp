#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE.
   */
   
   
  VectorXd rmse(4);
  rmse << 0,0,0,0;


  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  
  if ((estimations.size()==0) || (estimations.size()!=ground_truth.size())){
    std::cout<<"invalid inputs";
    return(rmse);
  }
    
  // TODO: accumulate squared residuals
  //cout<<estimations[0]<<endl;
  VectorXd a,b;
  VectorXd c;
  for (int i=0; i < estimations.size(); ++i) {
    
    a=estimations[i];
    b=ground_truth[i];
    c=a-b;
    c=(c.array()*c.array());
    rmse+=c;
    
  }

  // TODO: calculate the mean
    rmse=rmse/estimations.size();

  // TODO: calculate the squared root
  rmse=rmse.array().sqrt();

  // return the result
  return rmse;
}