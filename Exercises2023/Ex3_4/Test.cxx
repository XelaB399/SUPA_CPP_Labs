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
    
// Create a default FiniteFunction with no additions/changes to the default code
    FiniteFunction defaultFunction;
    defaultFunction.printInfo();
    defaultFunction.setRangeMax(4);
    defaultFunction.setRangeMin(-7);
    defaultFunction.plotFunction();
    defaultFunction.plotData(testData, 100, true);

// Create a NormalFunction with no additions/changes to the default code
    NormalFunction NormalFunction;
    NormalFunction.printInfo();
    NormalFunction.setRangeMax(4);
    NormalFunction.setRangeMin(-7); 
    NormalFunction.setMean(testData);
    NormalFunction.setStdDev(testData, (NormalFunction.getMean())); 
    NormalFunction.plotFunction();
    NormalFunction.plotData(testData, 100, true);

    return(0);
}