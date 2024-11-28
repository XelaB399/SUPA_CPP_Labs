//Alexandra Berger
//10.11.2024

#pragma once

//functions
std::vector<std::pair<double, double>> filereading(std::string filename);
void PrintVector(std::vector<std::pair<double, double>> data, int numberoflines);
void PrintVector(std::vector<double> data, int numberoflines);
std::vector<double> Magnitude(std::vector<std::pair<double, double>> data);
std::pair<double, double> FitLine(const std::vector<std::pair<double, double>> data, std::vector<std::pair<double, double>> Error);
std::vector<double> x_expy(const std::vector<std::pair<double, double>> data, int loopnum, std::vector<double> result);
void SaveResult(std::vector<double> data, std::string FunctionName);
