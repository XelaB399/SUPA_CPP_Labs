//Alexandra Berger
//10.11.2024

#include <iostream> //input output 
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include "CustomFunctions.h"

std::vector<std::pair<double, double>> filereading(std::string filename) {
// open a file in read mode.
   std::ifstream infile; 
   infile.open(filename); 
    //std::cout << "Reading from the file..." << std::endl;
//Declaring vector datastd which stores the x and y values
    std::vector<std::pair<double, double>> data;
//Declares sttring variable line to store individual lines read from file
    std::string line;
    //std::cout << "Placing x and y coords into a vector..." << std::endl;
// Read the file line by line
    while (std::getline(infile, line)) {
            // Output the line
            //std::cout << line << std::endl;
        std::istringstream iss(line);

        //Initialising x and y 
        double x, y; 
        char comma;
        // Extract x and y from the line, assuming they are separated by a comma
        if (iss >> x >> comma >> y && comma == ',') {
            data.emplace_back(x, y);
        }
    }
    // close the opened file.
   infile.close();
return data;
}

void PrintVector(std::vector<std::pair<double, double>> data, int numberoflines)
{
if (data.size() < numberoflines) {
    std::cout << "The number of lines in the file is less than the number of lines requested" << std::endl;
    std::cout << "Only first 5 lines will be printed" << std::endl;
    for (int i = 0; i < 5; i++) {
        std::cout <<"x, " << data[i].first << "   y, " << data[i].second << std::endl;
    }
}

else if (numberoflines == 0) {
    std::cout << "Showing all lines:" << std::endl;
    for (int i = 0; i < data.size(); i++) {
        std::cout <<"x, " << data[i].first << "   y, " << data[i].second << std::endl;
    }
}

else{
    std::cout << "Showing " << numberoflines << " lines:" << std::endl;
    for (int i = 0; i < numberoflines; i++) {
        std::cout <<"x, " << data[i].first << "   y, " << data[i].second << std::endl;
    }
}
return;
}

void PrintVector(std::vector<double> data, int numberoflines)
{
if (data.size() < numberoflines) {
    std::cout << "The number of lines in the file is less than the number of lines requested" << std::endl;
    std::cout << "Only first 5 lines will be printed" << std::endl;
    for (int i = 0; i < 5; i++) {
        std::cout <<"x^y, " << data[i] << std::endl;
    }
}
else if (numberoflines == 0) {
    std::cout << "Showing all lines:" << std::endl;
    for (int i = 0; i < data.size(); i++) {
        std::cout <<"x^y, " << data[i] << std::endl;
    }
}

else{
    std::cout << "Showing " << numberoflines << " lines:" << std::endl;
    for (int i = 0; i < numberoflines; i++) {
        std::cout <<"x^y, " << data[i] << std::endl;
    }
}
return;
}

std::vector<double> Magnitude(std::vector<std::pair<double, double>> data)
{
    std::vector<double> magnitudes;
    double mag;
    for (int i = 0; i < data.size(); i++) {
        mag = sqrt(pow(data[i].first, 2) + pow(data[i].second, 2));
        magnitudes.push_back(mag);
        //std::cout << "Magnitude of vectors, " << data[i].first << " and y, " << data[i].second << " is " << mag << std::endl;
    }
    return magnitudes;
}

std::pair<double, double> FitLine(const std::vector<std::pair<double, double>> data, std::vector<std::pair<double, double>> Error) {
    //Initialising variables
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    //num of data points
    double n = data.size();

    //Looping over each x and y pair in the data
    for (auto& pair : data){
        sum_x += pair.first;
        sum_y += pair.second;
        sum_x2 += pow(pair.first, 2);
        sum_xy += pair.first * pair.second;
    }

    //Performing least squares line fit where m is slope and c is y intercept
    double m = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    double c = (sum_y - m * sum_x) / n;


//Chi2 fit
 // Calculating expected values
    std::vector<std::pair<double, double>> expected;
    for (int i = 0; i < data.size(); i++) {
        double x = data[i].first;
        double y = m * x + c;
        expected.emplace_back(x, y);
    }

    // Initializing variables
    double chi2 = 0;

    // Looping over each y pair in the data
    for (int i = 0; i < data.size(); i++) {
        // Calculating the chi2 value for each data point
        double observed_y = data[i].second;
        double expected_y = expected[i].second;
        double error_y = Error[i].second;
        chi2 += std::pow((observed_y - expected_y), 2) / std::pow(error_y, 2);
    }

    //Creating results string
    std::ostringstream result_fl, result_chi2;
    result_fl << "Fitted line: y = " << m << "x + " << c;
    result_chi2 << "Chi2 value: " << chi2;

    //Printing results to terminal
    //std::cout << result_fl.str() << std::endl;
    //std::cout << result_chi2.str() << std::endl;

    //Saving results to a new file
    //Creating a new file to write the results to
    std::ofstream outfile("fitted_line_results.txt");
    //If file was created properly, write the results to the file
    if (outfile.is_open()) {
        outfile << result_fl.str() << std::endl;
        outfile << result_chi2.str() << std::endl;
        outfile.close();
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }

        return std::make_pair(m, c);
}


std::vector<double> x_expy(const std::vector<std::pair<double, double>> data, int loopnum, std::vector<double> result){
//std::vector<double> result; //Initialising result container
if (loopnum == data.size()) { //If all data has been analysed, return result and end
    //std::cerr << "Loop finished" << std::endl;
    return result;}

else{
    double x = data[loopnum].first;
    double y = data[loopnum].second;
    //std::cout << x << std::endl;
    //std::cout << y << std::endl;
    int y_round = std::round(y);
    //std::cout << y_round << std::endl;
    double xecpy = std::exp(y_round* std::log(x)); 
   // std::cout << xecpy << std::endl;
    result.push_back(xecpy);
    //std::cerr << result[loopnum] << std::endl;// Used for debugging


    return x_expy(data, loopnum + 1, result);
}
}


void SaveResult(std::vector<double> data, std::string FunctionName){
//Creating results string by looping through each data value
    std::ostringstream result;
    result << "Results for " << FunctionName << std::endl;
    for (int i = 0; i < data.size(); i++) {
      result<< data[i]<<"\n";
    }
    //Printing results to terminal for debugging
    //std::cout << result.str() << std::endl;

    //Creating a new file to write the results to
    std::ofstream outfile(FunctionName+".txt");
    //If file was created properly, write the results to the file
    if (outfile.is_open()) {
        outfile << result.str() << std::endl;
        std::cerr << "Results file created" << std::endl;
        outfile.close();
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }
    return;
}

//Fitted line will already save itself to a file when run, I have overloaded the 
//the print function line instead to show I know how to do it 

//Not all fringe cases have been fixed: For example not all required numerical user inputs are tested
// to make sure they are actually a number. Some cases are tested for proof of ability to, but time limiattions
// mean this was not implemented for all numerical inputs.

//Lack of git history as I am a silly goose and kept forgetting to commit and push.



