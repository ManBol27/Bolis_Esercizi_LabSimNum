#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../GenRand/random.h"

using namespace std;

/*******************S_N_printer*********************/
/*
Scrivo una funzione che presi M numeri casuali, costruisce K variabili casuali S_N = (r_1 + r_2 + ... + r_N)/N 
Per questo esercizio si considera N=1,2,10,100. Voglio dunque stampare un file con 4 colonne (una per ogni N)
ciascuna delle quali contenente K diverse estrazioni per S_N
*/

// A questa funzione viene passato un array di M numeri casuali estratti da una certa distribuzione
// e un intero K che rappresenta il numero di variabili casuali S_N da costruire
void S_N_printer(int K, const std::string& File_name, const std::vector<double>& r) {
   

   std::vector<int> valori_N = {1, 2, 10, 100};  // Definisco per quali N voglio costruire le variabili S_N

   std::vector<double> S_N1(K, 0), S_N2(K,0), S_N10(K,0), S_N100(K,0);  //vettori che conterranno le variabili casuali S_N
   
   for (int N : valori_N) {
      for(int i = 0; i < K; ++i) {
         double sum = 0.0;
         for(int j = i*N; j < i*N+N; ++j) {
            sum += r[j];
         }
         double S_N = sum / N;
         if (N == 1) S_N1[i] = S_N;
         else if (N == 2) S_N2[i] = S_N;
         else if (N == 10) S_N10[i] = S_N;
         else if (N == 100) S_N100[i] = S_N;
      }
   }

//Apro il file in cui voglio scrivere i risultati
   std::ofstream file(File_name);

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }

   file << setw(15) << "N=1" << setw(15) << "N=2" << setw(15) << "N=10" << setw(15) << "N=100" << "\n";
   for(int i=0; i<K; ++i) {
      file << setw(15) << S_N1[i] << setw(15) << S_N2[i] << setw(15) << S_N10[i] << setw(15) << S_N100[i] << "\n";
   }
//Chiudo il file
   file.close();
   }

/**************************************/

/*************************************/

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
  andrò a generare un numero M di numeri pseudo-casuali */
   const int M = 1000000; // Quanti numeri casuali devo generare
   const int K = 1E4;    // Numero di variabili casuali S_N da costruire

   //Genero un numero M di numeri pseudo-casuali da distribuzione uniforme
   std::vector<double> r_unif(M);
   for(int i = 0; i < M; ++i) {
      r_unif[i] = rnd.Rannyu();      //random generator uniforme
   }

   //Genero un numero M di numeri pseudo-casuali da distribuzione esponenziale
   std::vector<double> r_expo(M);
   for(int i = 0; i < M; ++i) {
      r_expo[i] = rnd.Expo(1);      //random generator esponenziale con lambda=1
   }

   //Genero un numero M di numeri pseudo-casuali da distribuzione lorentziana
   std::vector<double> r_lorentz(M);   
   for(int i = 0; i < M; ++i) {
      r_lorentz[i] = rnd.Lorentz(0, 1);      //random generator lorentziana con mu=0 e gamma=1
   }

   //Stampo file con S_N per N=1,2,10,100 con distrubuzione uniforme
   std::string File_name = "S_N_unif.txt";
   S_N_printer(K, File_name, r_unif);

   //Stampo file con S_N per N=1,2,10,100 con distrubuzione esponenziale
   File_name = "S_N_expo.txt";
   S_N_printer(K, File_name, r_expo);

   //Stampo file con S_N per N=1,2,10,100 con distrubuzione lorentziana
   File_name = "S_N_lorentz.txt";
   S_N_printer(K, File_name, r_lorentz);

   
   rnd.SaveSeed();
   return 0;
}