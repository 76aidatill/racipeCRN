#include "header.h"

using namespace Rcpp;

unsigned u_seed = 123;// std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 u_generator (u_seed);

//uniformly distributed random number generator in (0,1) range
std::uniform_real_distribution<double> u_distribution(0.0,1.0);

unsigned g_seed = 123456; //
std::mt19937_64 g_generator (g_seed);
std::normal_distribution<double> g_distribution(0.0,1.0);

//This function calculates the squared Euclidean distance between two vectors
double sum_delta (std::vector<double> &exprx1, 
                  std::vector<double> &exprx2, int numberReactions)
{
    double ssq = 0.0;
    
    for (int i=0;i<numberReactions;i++){
        ssq+=pow(exprx1[i]-exprx2[i],2);
    }
    return ssq;
}

//This function calculates the squared Euclidean norm of a vector.
double cal_norm(
    std::vector<double> &vector,
    const int &numberReactions
){
    double ssq = 0.0;
    for(int i=0;i<numberReactions;i++){
        ssq += pow(vector[i], 2);
    }
    return ssq;
}
