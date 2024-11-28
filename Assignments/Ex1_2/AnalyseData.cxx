//Alexandra Berger
//10.11.2024

#include <iostream> //input output 
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include "CustomFunctions.h"


int main() {

//Read file and store x and y values in a vector
std::vector<std::pair<double, double>> data =filereading("data/input2D_float.txt"); 
std::vector<std::pair<double, double>> error =filereading("data/error2D_float.txt"); 

std::cout <<"Hello there!"<< std::endl;
std::cout <<"For your convenience the relevant data and error files have already been loaded in :D" << std::endl;
bool endloop = false;
while (endloop == false){
   std::cout <<"Data analysis options:" << std::endl;
   std::cout <<"0) Print loaded data to the terminal" << std::endl;
   std::cout <<"1) Calculate the magnitude of each data point" << std::endl;
   std::cout <<"2) Fit a straight line to the data and calculate associated chi^2" << std::endl;
   std::cout <<"3) Calculate x^y for each data point" << std::endl;

   std::cout <<"Type the associated number of the desired action " << std::endl;
   int ActionNum;
   std::cin >> ActionNum;
   
   while (std::cin.fail()){
         std::cin.clear(); // Clear the error flag
         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore invalid input
         std::cout << "Invalid input. Please enter an integer: ";
         std::cout<<" " << std::endl;
         std::cin >> ActionNum;
   }


switch(ActionNum){
      case 0:
      {
         //Print data to terminal
         std::string PrintOrNo;
         std::cout <<"Would you like to print the data to the terminal? y/n"<<std::endl;
         std::cin >> PrintOrNo;
         std::cout <<" " << std::endl;
         if (PrintOrNo == "y"){
            int numlines; //Sees how many lines to print
            std::cout <<"How many lines do you want to print? (type 0 to show all lines)" << std::endl;
            std::cin >> numlines;
            std::cout <<" " << std::endl;

         if (std::cin.fail()) { //Make sure input is an integer
             std::cin.clear(); // Clear the error flag
             std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore invalid input
             std::cout << "Invalid input. Please enter an integer: ";
             break;
         } else {
               std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore any extra input
         }
            //Prints the data to terminal
            PrintVector(data, numlines);
         }
         else{
            std::cout <<"No information will be printed to terminal" << std::endl;
         }
         break;
      }

      case 1:
      {
         //Calculate magnitude of each data point
         std::cout <<"Saving Magnitude data to a file" << std::endl;
         std::vector<double> magnitude = Magnitude(data);
         SaveResult(Magnitude(data), "Magnitude");
         break;
      }
   
      case 2:
      {
         std::cout <<"Fitting the data" << std::endl;
         std::cout <<"A file will be created containing the fit values and the chi^2 value" << std::endl;
         FitLine(data, error);
         break;
      }
   
      case 3:
      {
         std::vector<double> xexpy = x_expy(data, 0, {});
         //Print result to terminal?
         std::string PrintOrNo;
         std::cout <<"Would you like to print the result to the terminal? y/n"<<std::endl;
         std::cin >> PrintOrNo;
         std::cout <<" " << std::endl;
         if (PrintOrNo == "yes" or PrintOrNo == "y"){
            int numlines;
            std::cout <<"How many lines do you want to print? (type 0 to shown all lines)" << std::endl;
            std::cin >> numlines;
            std::cout <<" " << std::endl;
            //Prints the result to terminal
            PrintVector(xexpy, numlines);
         }
         else{
            std::cout <<"No information will be printed to terminal" << std::endl;
         }
         std::cout <<" " << std::endl;
         std::cout <<"Saving x^y data to a file" << std::endl;
         SaveResult(xexpy, "x^y");
         break;
      }
   
      default:
      {
         std::cout <<"Invalid input, please try again" << std::endl;
         break;
      }
   } //End of Switch 

    std::string end;
    std::cout <<"Would you like to end the loop? yes/no" << std::endl;
    std::cin >> end;
    if (end == "yes" or end =="y"){
        endloop = true;
    }
    else{
        std::cout <<"Lets continue then!" << std::endl;
    }
} //End of while loop


return 0;
}

//Note: Read comments at the end of the CustomFunctions.cxx file for more information on the code
//for copy paste convience:
   // g++ AnalyseData.cxx CustomFunctions.cxx -o AnalyseData
   // ./AnalyseData