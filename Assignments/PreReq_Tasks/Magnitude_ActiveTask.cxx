#include <iostream>
#include <math.h>

int main(int argc, char** argv) {

double x = std::atof(argv[1]);
double y = std::atof(argv[2]);

double Magnitude = std::pow((std::pow(y,2) + std::pow(x,2)),0.5);
std::cout<< Magnitude;


return 0;
}