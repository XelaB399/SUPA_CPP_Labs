#include <iostream>
#include <math.h>

int main() {

 double x = 2.3;
 double y = 4.5;

 double Magnitude = std::pow((std::pow(y,2) + std::pow(x,2)), 0.5);

std::cout<< Magnitude;

return 0;
}