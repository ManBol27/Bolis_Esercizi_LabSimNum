//////////////////////////////////////// ESERCIZITAZIONE 2 ////////////////////////////////////////
// 1. Integrale Monte Carlo 1D con sampling uniforme
// 2. Integrale Monte Carlo 1D con importance sampling

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../GenRand/random.h"

using namespace std;

double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {
    if (n == 0)
	return 0.0;
    else
	return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

/****************************************************/


/******************FUNZIONE DA INTEGRARE - UNIF SAMPLING**********/

double f_integrand(double x) {

   double f =  M_PI/2 * cos(x*M_PI/2);
   return f;
}


/*************************************/

/******************FUNZIONE DA INTEGRARE - IMPORTANCE SAMPLING**********/

double f_integrand2(double x) {

   double p = 1 - (M_PI/2) * (x - 0.5);   
   double f =  M_PI/2 * cos(x*M_PI/2);;
   double g = f/p;
   return g;
}


/*************************************/

/******************** FUNZIONE AVE PRINTER **********/

/* Funzione per la medie cumulative dell'integrale con relativo errore - CASO UNIFORM SAMPLING*/

void Ave_printer(int N, int M, int L, const std::string& File_name, const std::vector<double>& r) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle      medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += f_integrand(r[k]);   
   }
     ave[i] = sum1/L;    //media della funzione f(x) sul blocco i-esimo
     av2[i] = pow(ave[i], 2);      //media della funzione f(x) al quadrato sul blocco i-esimo
       }

   for(int i = 0; i < N; ++i) {
     for(int j = 0; j <= i; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
   }
     sum_prog[i] /= (i + 1);      //media della f(x) sui primi i blocchi
     su2_prog[i] /= (i + 1);      //media della f(x)^2 sui primi i blocchi
     err_prog[i] = error(sum_prog, su2_prog, i);      //dev standard cumulativa sui primi i blocchi
    }

//Apro il file in cui voglio scrivere i risultati
   std::ofstream file(File_name);

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }

   file << setw(15) << "Media" << setw(15) << "Errore" << "\n";
   for(int i=0; i<N; ++i) {
      file << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << "\n";
   }
//Chiudo il file
   file.close();
   }

/**************************************/


/******************** FUNZIONE AVE PRINTER2 **********/

/* Funzione per la medie cumulative dell'integrale con relativo errore - CASO IMPORTANCE SAMPLING*/

void Ave_printer2(int N, int M, int L, const std::string& File_name, const std::vector<double>& r) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle      medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += f_integrand2(r[k]);   
   }
     ave[i] = sum1/L;    //media della funzione f(x) sul blocco i-esimo
     av2[i] = pow(ave[i], 2);      //media della funzione f(x) al quadrato sul blocco i-esimo
       }

   for(int i = 0; i < N; ++i) {
     for(int j = 0; j <= i; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
   }
     sum_prog[i] /= (i + 1);      //media della f(x) sui primi i blocchi
     su2_prog[i] /= (i + 1);      //media della f(x)^2 sui primi i blocchi
     err_prog[i] = error(sum_prog, su2_prog, i);      //dev standard cumulativa sui primi i blocchi
    }

//Apro il file in cui voglio scrivere i risultati
   std::ofstream file(File_name);

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }

   file << setw(15) << "Media" << setw(15) << "Errore" << "\n";
   for(int i=0; i<N; ++i) {
      file << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << "\n";
   }
//Chiudo il file
   file.close();
   }

/**************************************/

/*************** MAIN ****************/

int main (int argc, char *argv[]){

// Inizializzazione generatore numeri casuali
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../GenRand/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../GenRand/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
	 input >> property;
	 if( property == "RANDOMSEED" ){
	    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	    rnd.SetRandom(seed,p1,p2);
	 }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  /* Inizializzo variabili del mio problema: 
  andrò a generare un numero M di numeri casuali raggruppati in N blocchi  */
   const int M = 1000000; // Quanti numeri casuali devo generare
   const int N = 100;    // Numero di blocchi in cui li divido
   const int L = M / N;  // Numero di elementi in ogni blocco

   std::vector<double> r(M);
   for(int i = 0; i < M; ++i) {
      r[i] = rnd.Rannyu();       // genero M numeri casuali estratti uniformemente tra 0 e 1
   }


   std::ofstream out("random_numbers.txt"); //stampo file in cui compaiono tutti i valori estratti dalla particolare p(x)
   std::vector<double> r_2(M);
   for(int i = 0; i < M; ++i) {
      r_2[i] = rnd.MC_sampling();       // genero M numeri casuali estratti da una distribuzione di probabilità p(x) = 1 - (pi/2)(x - 1/2)
      out << r_2[i] << endl;
   }

   out.close();

   //Calcolo la media di f(x) (aka l'integrale richiesto) cumulativa con relativo errore usando uniform sampling
   string File_name = "MC_unif.txt";
   Ave_printer(N, M, L, File_name, r);
   
   //Calcolo la media di g(x) (aka l'integrale richiesto) cumulativa con relativo errore usando importance sampling
   File_name = "MC_importance.txt";
   Ave_printer2(N, M, L, File_name, r_2);
    
   
   rnd.SaveSeed();
   return 0;
}