#include "header.h"
#include <Rcpp.h>
#include <utility>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;


// [[Rcpp::export]]

int simulateCRN(Rcpp::IntegerMatrix stoichMatrix,
                Rcpp::IntegerMatrix reactantMatrix,
                Rcpp::List config, String outFileSC, String outFileParams,
                String outFileIC, String outFileConverge,
                const int stepper = 1){
    
        unsigned int seed =  static_cast<unsigned int>
        (Rcpp::sample(32000,1,true)(0));
        std::mt19937_64 g_generator (seed);
    
    // Initialize the network
    size_t numberReactions = stoichMatrix.ncol();
    size_t numberSpecies = stoichMatrix.nrow();

    NumericVector simulationParameters =
        as<NumericVector>(config["simParams"]);
    NumericVector hyperParameters =
        as<NumericVector>(config["hyperParams"]);
    LogicalVector options = as<LogicalVector>(config["options"]);

    size_t numModels = static_cast<size_t>(simulationParameters[0]);
    size_t numIC = static_cast<size_t> (simulationParameters[1]);
    size_t outputPrecision = static_cast<size_t> (simulationParameters[2]);

    double kMin = hyperParameters[0];
    double kMax = hyperParameters[1];
    double icMin = hyperParameters[2];
    double icMax = hyperParameters[3];

    bool integrate = options[0];
    bool genParams = options[1];
    bool genIC = options[2];

    std::string fileNameSC = outFileSC;
    std::string fileNameParam = outFileParams;
    std::string fileNameIC = outFileIC;
    std::string fileNameConverge = outFileConverge;

    //Create output files if not there already
    std::ofstream outSC(fileNameSC, std::ios::out);
    if(!outSC.is_open()) {     Rcout << "Cannot open output file.\n";
      return 1;}
    
    std::ofstream outConv(fileNameConverge, std::ios::out);
    if(!outConv.is_open()) {     Rcout << "Cannot open output file.\n";
      return 1;}
    
    std::ifstream inParams;
    std::ofstream outParam;

    if(genParams){
      outParam.open(fileNameParam,std::ios::out);
    }
    else {
      inParams.open(fileNameParam,
                         std::ifstream::in);
      if(!inParams.is_open()) {     Rcout <<fileNameParam
        << "Cannot open input file for reading parameters.\n";  return 1;
      }
    }

    std::ifstream inIC;
    std::fstream outIC;
    if(genIC){
      outIC.open(fileNameIC,std::ios::out);
    }
    else
    {
      inIC.open(fileNameIC,std::ifstream::in);
      if(!inIC.is_open()) {
        Rcout <<"Cannot open input file for reading initial conditions.\n";
        return 1;
        }
    }

        for(size_t modelCount=0;modelCount<numModels;modelCount++)
      {
            if((static_cast<size_t>(10*modelCount) % numModels) == 0){
              Rcout<<"====";
            }

            // Check for user interrupt and exit if there is any.
            if (modelCount % 100 == 0)
                Rcpp::checkUserInterrupt();

            //Initialize rate constants for reactions 
            std::vector<double> kConst(numberReactions);

            if(inParams.is_open())
                {
                for(size_t reacCount1=0;reacCount1<numberReactions;reacCount1++)
                    {    inParams >>  kConst[reacCount1];}
            }
            else
                {
                for(size_t reacCount1=0;reacCount1<numberReactions;reacCount1++){
                    kConst[reacCount1]=kMin+(kMax-kMin)*u_distribution(u_generator);
                }
                for(size_t reacCount1=0;reacCount1<numberReactions;reacCount1++){
                    outParam<<std::setprecision(outputPrecision)<<kConst[reacCount1]<<"\t";
                }
                outParam<<"\n";

            }

            for(size_t icCount=0;icCount<numIC;icCount++){
                std::vector <double> concSpecies(numberSpecies);
                //array for current species concentration
                std::vector <double> concSpecies0(numberSpecies);
                //array for initial species concentration

                if(!genIC)
                    {
                    for(size_t icCounter=0;icCounter <numberSpecies; icCounter++)
                        {
                        inIC >> concSpecies0[icCounter];
                    }
                }
                else {
                    for(size_t speciesCount=0; speciesCount<numberSpecies; speciesCount++){
                        concSpecies0[speciesCount] = icMin+(icMax-icMin)*u_distribution(u_generator);
                    }

                    for(size_t speciesCount=0; speciesCount<numberSpecies; speciesCount++){
                        outIC<<std::setprecision(outputPrecision)
                        <<concSpecies0[speciesCount]<<"\t";
                    }
                    outIC<<"\n";
                }

                for(size_t speciesCount=0; speciesCount<numberSpecies; speciesCount++)
                    {
                        concSpecies[speciesCount]=concSpecies0[speciesCount];
                }

                if(integrate){
                    //ADD INTEGRATION LATER
                    Rcout << "Integration currently unavailable";
                }

                outSC<<"\n";
            }
      }
    outSC.close();
    outParam.close();
    outIC.close();
    outConv.close();
    if(inIC.is_open()) inIC.close();
    if(inParams.is_open()) inParams.close();
    Rcout<<"\n";

    return 0;
}