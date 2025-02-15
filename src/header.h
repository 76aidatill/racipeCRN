#ifndef HEADER_H
#define HEADER_H

#include <Rcpp.h>
#include <iostream>
#include <random>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

//uniformly distributed random number generator in (0,1) range
// Shifted hill function
extern std::mt19937_64 u_generator;
extern std::uniform_real_distribution<double> u_distribution;
extern std::mt19937_64 g_generator;
// Gaussian distributed random number generator with mean 0 and 1
//standard deviation
extern std::normal_distribution<double> g_distribution;

double sum_delta (std::vector<double> &exprx1, 
                  std::vector<double> &exprx2, int numberReactions);

double cal_norm(std::vector<double> &vector, 
                const int &numberReactions);


#endif