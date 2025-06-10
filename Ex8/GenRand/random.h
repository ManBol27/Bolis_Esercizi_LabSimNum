/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Default constructor
  Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to save the seed to a file
  void SaveSeed();
  // Method to generate a random number in the range [0,1)
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);
  // Method to generate a random number with an exponential distribution
  double Expo(double lambda);
  // Method to generate a random number with a Lorentzian distribution  
  double Lorentz(double mu, double gamma);
  //Method used in Excercise 2
  double MC_sampling(void);
  //Method used in Excercise 2 for extracing x from p(x) = sin(x)
  double sinusoidal(void);
  //Metropolis sampling function with uniform transition
  std::vector<std::vector<double>> Metropolis_sampling (double (*p)(double, double, double, double, double), double S, int M, int extractions, std::string process, double mu, double sigma);
  //Metropolis sampling with Gauss transition
  std::vector<std::vector<double>> Metropolis_sampling_gauss (double (*p)(double, double, double), double S, int M, int extractions, std::string process);

};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
