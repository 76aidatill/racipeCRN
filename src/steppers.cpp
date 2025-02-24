#include "header.h"
#include <Rcpp.h>
using namespace Rcpp;

// This function evaluates the reaction rate vector at a given point in time
void calcReactVec(const size_t &numberReactions,
                  const size_t &numberSpecies,
                  std::vector<double> &speciesConc,
                  std::vector<double> &reactVector,
                  std::vector<double> &kConst,
                  Rcpp::IntegerMatrix reactantMatrix){
    double speciesFactor;
    for(size_t i=0; i<numberReactions; i++){
        reactVector[i] = kConst[i];
        for(size_t j=0; j<numberSpecies; j++){
            speciesFactor = pow(speciesConc[j], reactantMatrix( j , i ));
            reactVector[i] = reactVector[i]*speciesFactor;
        }
    }
}

// Calculate the derivative at a point in time.
void calcDeriv(const size_t &numberReactions,
                const size_t &numberSpecies,
                std::vector<double> &speciesConc,
                std::vector<double> &derivVec,
                std::vector<double> &reactVector,
                std::vector<double> &kConst,
                Rcpp::IntegerMatrix reactantMatrix,
                Rcpp::IntegerMatrix stoichMatrix){
    calcReactVec(numberReactions, numberSpecies, speciesConc, reactVector,
                kConst, reactantMatrix);
    
    for(size_t i=0; i<numberSpecies; i++){
        for(size_t j=0; j<numberReactions; j++){
            derivVec[i] += stoichMatrix( i , j ) * reactVector[j];
        }
    }
    }

// Use the Euler method to simulate an iteration, defined by some number of integration steps
void simIterE(const size_t &numberReactions,
            const size_t &numberSpecies,
            const double &h,
            const int &numStepsIteration,
            std::vector<double> &speciesConc,
            std::vector<double> &kConst,
            Rcpp::IntegerMatrix reactantMatrix,
            Rcpp::IntegerMatrix stoichMatrix){

    //Create derivative and reaction vector
    std::vector<double> derivVec(numberSpecies, 0.0);
    std::vector<double> reactVector(numberReactions, 0.0);

    for(int step=0; step<numStepsIteration; step++){
        calcDeriv(numberReactions, numberSpecies, speciesConc, derivVec,
            reactVector, kConst, reactantMatrix, stoichMatrix);
        
        for(size_t i=0; i<numberSpecies; i++){
            speciesConc[i] = speciesConc[i] + derivVec[i]*h;
            derivVec[i] = 0.0;
        }
    }
    
}

//Run iterations and report concentrations and convergence data
void runStepper(std::vector<double> &speciesConc,
        std::ofstream &outSC,
        std::ofstream &outConv,
        const size_t &outputPrecision,
        const size_t &numberReactions,
        const size_t &numberSpecies,
        const double &h,
        const double &convergThresh,
        const int &numStepsIteration,
        const int &numIterations,
        const int &stepper,
        std::vector<double> &kConst,
        Rcpp::IntegerMatrix reactantMatrix,
        Rcpp::IntegerMatrix stoichMatrix){

    // Make vector for comparisons between iterations
    std::vector<double> prevSpeciesConc(numberSpecies);

    bool isConverged = false;
    for(int testIter=0; testIter<numIterations; testIter++){
        for(size_t species=0; species<numberSpecies; species++){
            prevSpeciesConc[species] = speciesConc[species];
        }
        switch(stepper){
            case 1:
            // Euler Method
            simIterE(numberReactions, numberSpecies, h, numStepsIteration,
                    speciesConc, kConst, reactantMatrix, stoichMatrix);
            break;

            default:
            Rcout<< "Error in specifying the stepper.\n";
        }
        double test_delta;
        test_delta = sum_delta(prevSpeciesConc, speciesConc, numberSpecies);
        if(test_delta < convergThresh){
            for(size_t speciesCount=0; speciesCount<numberSpecies; speciesCount++){
                outSC<<std::setprecision(outputPrecision)
                <<speciesConc[speciesCount]<<"\t";
            }
            isConverged = true;
            outConv<<isConverged<<"\t" << testIter <<"\n";
            break;
        }
    }
    if(isConverged == false){
        for(size_t speciesCount=0; speciesCount<numberSpecies; speciesCount++){
            outSC<<std::setprecision(outputPrecision)
            <<speciesConc[speciesCount]<<"\t";
        }
        outConv<<isConverged<<"\t" << numIterations <<"\n";
    }
}