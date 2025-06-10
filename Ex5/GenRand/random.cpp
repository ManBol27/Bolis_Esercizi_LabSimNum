/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ostream>
#include <string>
#include "random.h"

using namespace std;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

double Random :: Expo (double lambda){
  // This function generates a random number from an exponential distribution with given lambda
  double U = Rannyu();
  return -log(1-U)/lambda;
}

double Random :: Lorentz (double mu, double gamma){
  // This function generates a random number from a Lorentzian distribution with given mu and gamma
  double U = Rannyu();
  return mu + gamma*tan(M_PI*(U-0.5));
}

//Distribuzione usata nell'esercizio 2 per Monte Carlo importance sampling
double Random :: MC_sampling (void){
  // This function generates a random number from a Lorentzian distribution with given mu and gamma
  double U = Rannyu();
  return ((4+M_PI) - sqrt(pow((4 +M_PI),2) - 16*M_PI*U))/(2*M_PI);
}

//Distribuzione usata nell'esercizio 2 per sampling da p(x) = sin(x)
double Random :: sinusoidal (void){
  // This function generates a random number from a Lorentzian distribution with given mu and gamma
  double U = Rannyu();
  return acos(1 - 2*U);
}


/**************************************************************************************/
/****************************FUNZIONE ALGORITMO METROPOLIS*****************************/
/**************************************************************************************/

//Funzione per fare un sampling di una qualsiasi funzione p(x,y,z) con algoritmo metropolis
//devo passare la funzione p(x,y,z) come argomento
//ritorna una posizione x,y,z
std::vector<std::vector<double>>Random :: Metropolis_sampling_gauss (double (*p)(double, double, double), double S, int M, int extractions, std::string process){
  
  std::vector<double> passo = {0.0, 0.0, 0.0};  //creo un vettore step con modulo S uniforme di 3 dimensioni 
  //Iniziallizzo il vettore che conterrà le estrazioni in 3d dopo la fase di equilibrazione
  std::vector<std::vector<double>> X_sampled(extractions, std::vector<double>(3, 0.0)); //il secondo indice ha 3 valori tutti iniziallizzati a 0.
  
  //IMPOSTO LA POSIZIONE INIZIALE
  std::vector<double> X_position = {0.0, 0.0, 0.0};  //vettore posizione, è la X_n che evolve passo dopo passo
 

  for (int i=0; i<M; i++){ //ciclo per M passi
    passo [0] = Gauss(0.,pow(S,2));
    passo [1] = Gauss(0.,pow(S,2));
    passo [2] = Gauss(0.,pow(S,2));  
    //scrivo la funzione di accettazione, ipotesi T(x,y) = T(y,x) (simmetria della transizione proposta)
    double A = std :: min(1., p(X_position[0]+passo[0],X_position[1]+passo[1],X_position[2]+passo[2])/p(X_position[0],X_position[1],X_position[2]));
    double r = Rannyu();
    if (r < A) { //accetto il passo
      X_position[0] += passo[0];
      X_position[1] += passo[1];
      X_position[2] += passo[2];
    }
    //Altrimenti non accetto il passo e rimango nella posizione vecchia
  }

  //Dopo la fase di equilibrazione, inizio a campionare la funzione p(x,y,z), l'ultimo valore della posizione è il primo che campiono
  X_sampled[0][0] = X_position[0];
  X_sampled[0][1] = X_position[1];
  X_sampled[0][2] = X_position[2];

  // int j=0; //contatore per quante estrazioni faccio, indipendentemente che vengano accettate o meno

  //Ora ciclo per le estrazioni vere e proprie (posizioni in 3d a partire dalla PDF passata come argomento)
  int accepted = 0;
  for (int i=1; i<extractions; i++){ //parto dall'indice i = 1 perchè la posizione 0 è gia settata nelle righe sopra

    passo [0] = Gauss(0.,pow(S,2));
    passo [1] = Gauss(0.,pow(S,2));
    passo [2] = Gauss(0.,pow(S,2)); 
    //scrivo la funzione di accettazione, ipotesi T(x,y) = T(y,x) (simmetria della transizione proposta)
    double A = std :: min(1., p(X_sampled[i-1][0]+passo[0],X_sampled[i-1][1]+passo[1],X_sampled[i-1][2]+passo[2])/p(X_sampled[i-1][0],X_sampled[i-1][1],X_sampled[i-1][2]));
    double r = Rannyu();
    if (r < A) { //accetto il passo
      X_sampled[i][0] = X_sampled[i-1][0] + passo[0];
      X_sampled[i][1] = X_sampled[i-1][1] + passo[1];
      X_sampled[i][2] = X_sampled[i-1][2] + passo[2];
      accepted++;
    }
    else { //non accetto il passo e rimango nella posizione vecchia
      X_sampled[i][0] = X_sampled[i-1][0];
      X_sampled[i][1] = X_sampled[i-1][1];
      X_sampled[i][2] = X_sampled[i-1][2];
    }
  }
  
  cout<<"Current acceptance rate: "<<double(accepted)/extractions <<endl;

  return X_sampled; //ritorno il vettore con le estrazioni in 3d
}


/**************************************************************************************/
/****************************END METROPOLIS********************************************/
/**************************************************************************************/

/**************************************************************************************/
/****************************FUNZIONE ALGORITMO METROPOLIS*****************************/
/**************************************************************************************/

//Funzione per fare un sampling di una qualsiasi funzione p(x,y,z) con algoritmo metropolis
//devo passare la funzione p(x,y,z) come argomento
//ritorna una posizione x,y,z
std::vector<std::vector<double>>Random :: Metropolis_sampling (double (*p)(double, double, double), double S, int M, int extractions, std::string process){
  
  std::vector<double> passo = {0.0, 0.0, 0.0};  //creo un vettore step con modulo S uniforme di 3 dimensioni 
  //Iniziallizzo il vettore che conterrà le estrazioni in 3d dopo la fase di equilibrazione
  std::vector<std::vector<double>> X_sampled(extractions, std::vector<double>(3, 0.0)); //il secondo indice ha 3 valori tutti iniziallizzati a 0.
  std::vector<double> X_position = {0., 0., 0.};  //vettore posizione, è la X_n che evolve passo dopo passo
 

  for (int i=0; i<M; i++){ //ciclo per M passi
    passo [0] = S * Rannyu(-1.,1.);
    passo [1] = S * Rannyu(-1.,1.);
    passo [2] = S * Rannyu(-1.,1.);  
    //scrivo la funzione di accettazione, ipotesi T(x,y) = T(y,x) (simmetria della transizione proposta)
    double A = std :: min(1., p(X_position[0]+passo[0],X_position[1]+passo[1],X_position[2]+passo[2])/p(X_position[0],X_position[1],X_position[2]));
    double r = Rannyu();
    if (r < A) { //accetto il passo
      X_position[0] += passo[0];
      X_position[1] += passo[1];
      X_position[2] += passo[2];
    }
    //Altrimenti non accetto il passo e rimango nella posizione vecchia
  }

  //Dopo la fase di equilibrazione, inizio a campionare la funzione p(x,y,z), l'ultimo valore della posizione è il primo che campiono
  X_sampled[0][0] = X_position[0];
  X_sampled[0][1] = X_position[1];
  X_sampled[0][2] = X_position[2];

  // int j=0; //contatore per quante estrazioni faccio, indipendentemente che vengano accettate o meno

  //Ora ciclo per le estrazioni vere e proprie (posizioni in 3d a partire dalla PDF passata come argomento)
  int accepted = 0;
  for (int i=1; i<extractions; i++){ //parto dall'indice i = 1 perchè la posizione 0 è gia settata nelle righe sopra

    passo [0] = S * Rannyu(-1.,1.);
    passo [1] = S * Rannyu(-1.,1.);
    passo [2] = S * Rannyu(-1.,1.);  
    //scrivo la funzione di accettazione, ipotesi T(x,y) = T(y,x) (simmetria della transizione proposta)
    double A = std :: min(1., p(X_sampled[i-1][0]+passo[0],X_sampled[i-1][1]+passo[1],X_sampled[i-1][2]+passo[2])/p(X_sampled[i-1][0],X_sampled[i-1][1],X_sampled[i-1][2]));
    double r = Rannyu();
    if (r < A) { //accetto il passo
      X_sampled[i][0] = X_sampled[i-1][0] + passo[0];
      X_sampled[i][1] = X_sampled[i-1][1] + passo[1];
      X_sampled[i][2] = X_sampled[i-1][2] + passo[2];
      accepted++;
    }
    else { //non accetto il passo e rimango nella posizione vecchia
      X_sampled[i][0] = X_sampled[i-1][0];
      X_sampled[i][1] = X_sampled[i-1][1];
      X_sampled[i][2] = X_sampled[i-1][2];
    }
  }
  
  cout<<"Current acceptance rate: "<<double(accepted)/extractions <<endl;

  return X_sampled; //ritorno il vettore con le estrazioni in 3d
}


/**************************************************************************************/
/****************************END METROPOLIS********************************************/
/**************************************************************************************/

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
