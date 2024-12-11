#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <random>
#include <algorithm> // for std::min

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}
NormalFunction::NormalFunction() : FiniteFunction() {
    m_stddev = 1.0;
    m_mean = 0.0;
    this->checkPath("NormalFunction");
}
CauchyLorentzFunction::CauchyLorentzFunction() : FiniteFunction() {
    m_x0 = 0.0;
    m_gamma = 1.0;
    this->checkPath("CauchyLorentzFunction");
}
NegativeCrystalBallFunction::NegativeCrystalBallFunction() : NormalFunction() {
    m_alpha = 1.0;
    m_n = 1.0;
    this->checkPath("NegativeCrystalBallFunction");
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

NormalFunction::NormalFunction(double min_range, double max_range, std::string outfile, double stddev, double mean){
  m_RMin = min_range;
  m_RMax = max_range;
  m_stddev = stddev;
  m_mean = mean;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

CauchyLorentzFunction::CauchyLorentzFunction(double range_min, double range_max, double gamma, double x0, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  m_x0 = x0;
  m_gamma = gamma;
  this->checkPath(outfile); //Use provided string to name output files
}
NegativeCrystalBallFunction::NegativeCrystalBallFunction(double min_range, double max_range, std::string outfile, double alpha, double n, double stddev, double mean){
  m_RMin = min_range;
  m_RMax = max_range;
  m_alpha = alpha;
  m_n = n;
  m_stddev = stddev;
  m_mean = mean;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}



//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}


/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};

void NormalFunction::setStdDev(std::vector<double> data, double mean) {
    double variance = 0.0;
    for (double value : data) {
        variance += std::pow(value - mean, 2);
    }
  m_stddev = std::sqrt(variance / data.size());}

void NormalFunction::setMean(std::vector<double> data) {
    double sum = 0.0;
    for (double value : data) {sum += value;}
  m_mean = sum / data.size();}

void CauchyLorentzFunction::setx0(double x0) {m_x0 = x0;};
void CauchyLorentzFunction::setGamma(double gamma) {m_gamma = gamma;};

void NegativeCrystalBallFunction::setAlpha(double alpha) {m_alpha = alpha;};
void NegativeCrystalBallFunction::setN(double n) {m_n = n;};

void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

double NormalFunction::getStdDev() {return m_stddev;};
double NormalFunction::getMean() {return m_mean;};

double CauchyLorentzFunction::getX0() {return m_x0;};
double CauchyLorentzFunction::getGamma() {return m_gamma;};

double NegativeCrystalBallFunction::getN() {return m_n;};
double NegativeCrystalBallFunction::getAlpha() {return m_alpha;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

//Normal Distribution
double NormalFunction::normalDistribution(double x) {
  //Getting variables
double m_stddev = this->getStdDev();
double m_mean = this->getMean();
  //Checking if the standard deviation is zero, return 0 if is
  if (m_stddev == 0) {return 0;}
  //Calculating the normal distribution in parts for convenience
double A = 1/(m_stddev*sqrt(2*M_PI));
double B = ((x-m_mean)*(x-m_mean))/(m_stddev*m_stddev);
double C = exp(-0.5*B);
return A*C;}
double NormalFunction::callFunction(double x) {return this->normalDistribution(x);}; 


//Cauchy-Lorentz Distribution
double CauchyLorentzFunction::cauchylorentzDistribution(double x) {
  //Getting variables
double m_x0 = this->getX0();
double m_gamma = this->getGamma();
  //Checking if the gamma is zero, return 0 if is
if (m_gamma == 0) {return 0;}
  //Calculating the Cauchy-Lorentz distribution 
double A = ((x-m_x0)*(x-m_x0))/m_gamma;
double B = (M_PI*m_gamma*(1+A));
  return 1/B;};
double CauchyLorentzFunction::callFunction(double x) {return this->cauchylorentzDistribution(x);}; //(overridable)

double NegativeCrystalBallFunction::NegativeCrystalBallDistribution(double x) {
  //Error checking
    if (m_stddev == 0.0) {
        // Make known if error
        throw std::invalid_argument("Standard deviation must be non-zero.");
    }
    if (m_n <= 1.0) {
        // Make known if error
        throw std::invalid_argument("n must be greater than one.");
    }
    if (m_alpha <= 0.0) {
        //Make known if error
        throw std::invalid_argument("Alpha must be greater than zero.");
    }

    //Calculating the negative crystal ball function in parts for convenience
    double A = pow(m_n / m_alpha, m_n) * exp(-0.5 * pow(m_alpha, 2));
    double B = m_n / m_alpha - m_alpha;
    double C = m_n / m_alpha * (1.0 / (m_n - 1)) * exp(-0.5 * pow(m_alpha, 2));
    double D = sqrt(M_PI / 2) * (1 + erf(m_alpha / sqrt(2)));
    double N = 1.0 / (m_stddev * (C + D));
    if ((x - m_mean)/m_stddev <= -m_alpha) {
        return N * A * pow(B - (x - m_mean) / m_stddev, -m_n);
    }
    else {
        return N * exp(-0.5 * pow((x - m_mean) / m_stddev, 2));
    }
}
double NegativeCrystalBallFunction::callFunction(double x) {return this->NegativeCrystalBallDistribution(x);}; //(overridable)



/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //trapezium rule 
  double step = (m_RMax - m_RMin)/(double)Ndiv; //Step size for integration (h)
  double x = m_RMin;                            //Start at the lower bound
  double sum = 0;                               //Initialising sum
  for (int i = 0; i < Ndiv; i++){               //Loop over the number of divisions
    sum += ((this->callFunction(x) + this->callFunction(x+step))/2)*step;   //Add the area of the trapezium to the total
    x = x + step;
  }
  return sum;                                //Return the total area under the curve
}

double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

void NormalFunction::printInfo(){
  std::cout << "StdDiv: " << m_stddev << std::endl;
  std::cout << "Mean: " << m_mean << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

void CauchyLorentzFunction::printInfo(){
  std::cout << "X0: " << m_x0 << std::endl;
  std::cout << "Gamma: " << m_gamma << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

void NegativeCrystalBallFunction::printInfo(){
  std::cout << "StdDiv: " << m_stddev << std::endl;
  std::cout << "Mean: " << m_mean << std::endl;
  std::cout << "n: " << m_n << std::endl;
  std::cout << "Alpha: " << m_alpha << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}




//Metropolis Algorithm//

void FiniteFunction::metropolisSampling(int numSamples, float posStd) {
    //Initialising the random number generator and associated variables
    std::vector<double> sampleData; 
    std::vector<double> randomX(numSamples); 
    std::vector<double> randomY(numSamples);
    std::vector<double> randomT(numSamples);
    int randomNum = 1;

    std::uniform_real_distribution<double> uniformDistT{0.0, 1.0}; // Uniform distribution for T
    std::random_device rd;
    std::mt19937 mtEngine{rd()}; // Mersenne Twister engine random number generator (found online) 
    std::uniform_real_distribution<double> rndNumber{m_RMin, m_RMax}; //range of the function uniform dist
//Generating the initial random sample    
    for (int j = 0; j < randomNum; j++) {
        // Step 1: Generate an initial random sample 
        double rndX = rndNumber(mtEngine);
        randomX[j] = rndX;
            //std::cout << "Initial Random X[0]: " << randomX[0] << std::endl;
    }
//Looping through the number of samples
    for (int i = 0; i < numSamples; i++) {
// Step 2: Generate a second random number 'y' from a normal distribution centered on 'xi'
        float centre = randomX[i];
        std::normal_distribution<float> normalPDF{centre, posStd};
        double rndY = normalPDF(mtEngine);
        randomY[i] = rndY;
       // std::cout << "Random Y[" << i << "]: " << rndY << std::endl;

//Step 3: Computing A = min(f(y)/f(x), 1)
        double f_x = callFunction(randomX[i]);
        double f_y = callFunction(randomY[i]); 
        double A = std::min(f_y / f_x, 1.0);
       // std::cout << "A[" << i << "]: " << A << std::endl;

//Step 4: Generate a random number 'T' between 0 and 1
        double T = uniformDistT(mtEngine);
        randomT[i] = T;
        // std::cout << "Random T[" << i << "]: " << T << std::endl;

// Step 4.5 and 5: If 'T < A', then accept 'y'.  If accepted y, set xi+1=y, otherwise xi+1=xi
        if (randomT[i] < A) { //Accept Y
          randomX[i + 1] = randomY[i];
           //  std::cout << "Accepted Y[" << i << "]: " << rndY << std::endl;
        } else { //Reject Y
          randomX[i + 1] = randomX[i];
          //  std::cout << "T[" << randomT[i] << "],  A [" << A << "]: " << rndX << std::endl;
        }
        //Storing data for plotting
        sampleData.push_back(randomX[i + 1]);
    }

    this->plotData(sampleData, 100, false);

    // Save the sampled data to a .txt file
    std::ofstream outFile("Metrosampled_data_Example.txt");
    if (outFile.is_open()) {
        for (const auto& value : sampleData) {
            outFile << value << "\n";
        }
        outFile.close();
        std::cout << "Sampled data saved to Metrosampled_data_Example.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }
  

    }






/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}
