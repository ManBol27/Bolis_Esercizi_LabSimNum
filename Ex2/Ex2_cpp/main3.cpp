//////////////////////////////////////// ESERCIZITAZIONE 2 ////////////////////////////////////////
// 2. Random Walk 3D continuo

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../GenRand/random.h"

using namespace std;


/****************************************************/


/******************** FUNZIONE AVE PRINTER **********/

/* Funzione per la medie cumulative di pigreco con relativo errore */

void RW_printer(int N, int B, int Z, const std::string& File_name, const std::vector<double>& moves) {
   
   std::vector<std::vector<double>> mean_matrix(N, std::vector<double>(B, 0.0));    /*matrice dove i-esima riga si riferisce alla distanza calcolata dopo i-passi, 
   mentre la s-esima colonna indica il blocco da 100 simulazioni su cui sto mediando*/
   std::vector<std::vector<double>> mean_matrix_squared(N, std::vector<double>(B, 0.0));     //idem per i quadrati 

   std::vector<double> sum_block(N), sum2_block(N), err(N); 
   
   //In questo caso N indica il numero di passi, mi serve considerare 1E4 simulazioni da N=1, 2, 3....., 100 passi ciasuna
   for(int i = 0; i < N; ++i) {  
     for (int s = 0; s < B; s++) {   //B numero di blocchi in cui divido le 10^4 simulazioni (B=100 blocchi da Z=100 simulazioni ciascuno)
         double sum1 = 0.0;   
         for(int j = 0; j < Z; ++j) {
            int k = i + s * N * N + j * N;   //struttura diversa rispetto agli esercizi precedenti!
            sum1 += moves[k];   
         }
         mean_matrix[i][s] = sum1/Z;        //media delle posizioni dopo un numero di i passi per il blocco s-esimo 
         mean_matrix_squared[i][s] = pow(mean_matrix[i][s], 2);      //matrice dei quadrati delle medie
     }

   }

   for(int i = 0; i < N; ++i) {
     for(int s = 0; s < B; ++s) {
      sum_block[i] += mean_matrix[i][s];
      sum2_block[i] += mean_matrix_squared[i][s];
   }
     sum_block[i] /= B;      //media sui B blocchi delle distanze raggiunte dopo numero i+1 di passi
     sum2_block[i] /= B;      //idem al quadrato
     err[i] = sqrt((sum2_block[i]-pow(sum_block[i], 2))/B);    //dev standard calcolata sui blocchi per la distanza raggiunta mediamente dai random walks dopo i+1 passi
    }

//Apro il file in cui voglio scrivere i risultati
   std::ofstream file(File_name);

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }


// In fase di stampa siccome devo plottare la radice della distanza quadratica media, stampo anche l'errore associato (calcolato con propagazione degli errori)
   file << setw(15) << "Media" << setw(15) << "Errore" << "\n";
   for(int i=0; i<N; ++i) {
      file << setw(15) << sqrt(sum_block[i]) << setw(15) << err[i]/(2*sqrt(sum_block[i])) << "\n";   
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
   const int M = 1000000; // Quanti numeri casuali devo generare tra 0 e 5
   const int N = 100;    // numero passi di un singolo random walk
   const int L = M / N;  // Numero di simulazioni (numero di random walks con 100 passi ciascuno)
   const int B = 100;     //Numero di blocchi in cui separo le 10000 simulazioni per l'i-esimo numero di passi
   const int Z = 100;      //Numero di simulazioni in ogni blocco per l'i-esimo numero di passi

   //Genero il vettore delle estrazioni casuali tra 0 e 5 
   //Considero le seguenti corrispondenze:
   // 0 -> +x 1 -> -x 2 -> +y 3 -> -y 4 -> +z 5 -> -z   
   std::vector<double> phi_angle(M);
   std::vector<double> theta_angle(M);
   for(int i = 0; i < M; ++i) {
      phi_angle[i] = rnd.Rannyu(0,2*M_PI);       // genero M valori per phi compresi tra 0 e 2pi
      theta_angle[i] = rnd.sinusoidal();
   }

   std::vector<double> moves(M);
   for(int j = 0; j < L; ++j) {
      double x_pos = 0;
      double y_pos = 0;
      double z_pos = 0;
      double distance_squared = 0;
      for(int i=0; i<N; ++i) {
         int k = i + j * N;
         x_pos += sin(theta_angle[k]) * cos(phi_angle[k]);
         y_pos += sin(theta_angle[k]) * sin(phi_angle[k]);
         z_pos += cos(theta_angle[k]);
         distance_squared = x_pos * x_pos + y_pos * y_pos + z_pos * z_pos;
         moves[k] = distance_squared;      // leggi la struttura del vettore moves dal notebook jupyter
      }
   }
      

   //Calcolo la media di f(x) cumulativa con relativo errore usando uniform sampling
   string File_name = "RW_continuo.txt";
   RW_printer(N, B, Z, File_name, moves);
   
   rnd.SaveSeed();
   return 0;
}