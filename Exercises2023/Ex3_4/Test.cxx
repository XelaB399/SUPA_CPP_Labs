#include <iostream>
#include "FiniteFunctions.h"
#include <fstream>
#include <sstream>
// #include <algorithm>
// #include <cmath>

int main() {

 // Read in the mystery data
    std::ifstream MysteryData;
    MysteryData.open("Outputs/data/MysteryData00110.txt");
//Making sure file opens correctly
    if (!MysteryData.is_open()) {
        std::cout << "Error opening file" << std::endl;
        exit(1);
    }
    std::vector<double> testData;
    std::string line;
    double dataLine;
//Read the data from the file
    while (std::getline(MysteryData, line)) {
        // Use an istringstream to parse the values separated by the given delimiter
        std::istringstream ss(line);
        std::string token;
        // Read the values separated by the specified delimiter
        std::getline(ss, token, ',');
        dataLine = std::stod(token); // Convert the string to the template type
        // Adding the values to the arrays
        testData.push_back(dataLine);
    }
    MysteryData.close();
    

// Metropolis Sampling
   FiniteFunction finiteFunction;
   finiteFunction.setRangeMin(-10);
   finiteFunction.setRangeMax(10);
  // finiteFunction.metropolisSampling(1000, 1.49717);

// Create a default FiniteFunction with no additions/changes to the default code
    FiniteFunction defaultFunction;
    defaultFunction.printInfo();
    defaultFunction.setRangeMax(14);
    defaultFunction.setRangeMin(-17);
    defaultFunction.plotFunction();
    defaultFunction.metropolisSampling(10000, 1.49717);
    defaultFunction.plotData(testData, 100, true);

// Create a NormalFunction with no additions/changes to the default code
    NormalFunction NormalFunction;
    NormalFunction.printInfo();
    NormalFunction.setRangeMax(14);
    NormalFunction.setRangeMin(-17); 
    NormalFunction.setMean(testData);
    NormalFunction.setStdDev(testData, (NormalFunction.getMean())); 
    NormalFunction.plotFunction();
    NormalFunction.metropolisSampling(10000, 1.49717);
    NormalFunction.plotData(testData, 100, true);

// Create a CauchyLorentzFunction with no additions/changes to the default code
    CauchyLorentzFunction CauchyLorentzFunction;
    CauchyLorentzFunction.setRangeMax(14);
    CauchyLorentzFunction.setRangeMin(-17);
    CauchyLorentzFunction.setGamma(1.6);
    CauchyLorentzFunction.setx0(-2);
    CauchyLorentzFunction.printInfo();
    CauchyLorentzFunction.plotFunction();
    CauchyLorentzFunction.metropolisSampling(10000, 1.49717);
    CauchyLorentzFunction.plotData(testData, 100, true);

// Create a NegativeCrystalBallFunction with no additions/changes to the default code
    NegativeCrystalBallFunction NegativeCrystalBallFunction;
    NegativeCrystalBallFunction.setRangeMax(14);
    NegativeCrystalBallFunction.setRangeMin(-17);
    NegativeCrystalBallFunction.setMean(testData);
    NegativeCrystalBallFunction.setStdDev(testData, (NegativeCrystalBallFunction.getMean()));
    NegativeCrystalBallFunction.setAlpha(2.4);
    NegativeCrystalBallFunction.setN(2);
    NegativeCrystalBallFunction.printInfo();
    NegativeCrystalBallFunction.plotFunction();
    NegativeCrystalBallFunction.plotData(testData, 100, true);
    NegativeCrystalBallFunction.metropolisSampling(10000, 1.49717);





    return(0);
}