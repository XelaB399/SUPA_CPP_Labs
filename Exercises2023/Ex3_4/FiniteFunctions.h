#include <string>
#include <vector>
#include <cmath>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

//////////////////////////////////////// 
////////////////////////////////////////
class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor

  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  

  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  void metropolisSampling(int numSamples, float proposalStd);

  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  
private:
  double invxsquared(double x); //The default functional form
};
//////////////////////////////////////// 
////////////////////////////////////////

//Defining Child classes of FiniteFunction
class NormalFunction : public FiniteFunction {
public:
    //Constructors
    NormalFunction(); // Default constructor
    NormalFunction(double min_range, double max_range, std::string outfile, double stddev, double mean); // Variable constructor 
    
    //Setters
    void setStdDev(std::vector<double> data, double mean);
    void setMean(std::vector<double> data);

    virtual double callFunction(double x) override; // overriding the callFunction from FiniteFunction
    virtual void printInfo() override; // overriding the printInfo from FiniteFunction

    //Getters
    double getStdDev();
    double getMean();

protected:
    double m_stddev;
    double m_mean;

private:
    double normalDistribution(double x);
};


class CauchyLorentzFunction : public FiniteFunction {
public:
    //Constructors
    CauchyLorentzFunction(); // Default constructor
    CauchyLorentzFunction(double min_range, double max_range, double gamma, double x0, std::string outfile); // Variable constructor 
    
    //Setters
    void setx0(double x0);
    void setGamma(double gamma);

    virtual double callFunction(double x) override; // overriding the callFunction from FiniteFunction
    virtual void printInfo() override; // overriding the printInfo from FiniteFunction

    //Getters
    double getX0();
    double getGamma();

protected:
    double m_x0;
    double m_gamma;

private:
    double cauchylorentzDistribution(double x);
};



class NegativeCrystalBallFunction : public NormalFunction {
public:
    NegativeCrystalBallFunction(); // Default constructor
    NegativeCrystalBallFunction(double min_range, double max_range, std::string outfile, double alpha, double n, double stddev, double mean); // Parameterized constructor
//   ~NegativeCrystalBallFunction(); //Destructor

    void setAlpha(double alpha);
    void setN(double n);

    virtual double callFunction(double x) override; 
    virtual void printInfo() override;

    double getAlpha();
    double getN();

protected:
    double m_alpha;
    double m_n;
  
private:
    double NegativeCrystalBallDistribution(double x);
};