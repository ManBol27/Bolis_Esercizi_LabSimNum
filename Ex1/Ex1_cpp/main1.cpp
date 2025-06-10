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
/******************** FUNZIONE AVE PRINTER **********/

/* Funzione per la medie cumulative con dev std cumulativa, per generico r-s
dove metterò s=0 per la media, e s=1/2 per la varianza */

void Ave_printer(int N, int M, int L, const std::string& File_name, double s, const std::vector<double>& r) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle      medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += pow(r[k]-s, 1+2*s);   // s=0 nel caso voglio la media, s=1/2 per la varianza
   }
     ave[i] = sum1 / L;      //media dell'i-esimo blocco
     av2[i] = pow(ave[i], 2);      //media al quadrato dell'i-esimo blocco
       }

   for(int i = 0; i < N; ++i) {
     for(int j = 0; j <= i; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
   }
     sum_prog[i] /= (i + 1);      //media delle medie sui primi i blocchi
     su2_prog[i] /= (i + 1);      //media delle medie quadrate sui primi i blocchi
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

/******************FUNZIONE CHI-QUADRO**********/

void Chi_quadro(int N, int L, int n_bins, const std::string& File_name, const std::vector<double>& r) {
   
   std::vector<double> chi2(N, 0);  //vettore che conterrà il chi-quadro per ogni blocco
   //Per ogni blocco, conto quanti elementi sono nell' n-esimo bin:
   for(int i = 0; i < N; ++i) {
      std::vector<int> V(n_bins, 0);   //vettore che conterrà il numero di elementi in ogni bin
      for(int j = 0; j < L; ++j) {
       int k = j + i * L;
       V[trunc(r[k]*100)]++;   /*se divido l'intervallo [0,1] in 100 bins, ogni bin avrà ampiezza 0.01. 
       Quindi moltiplico per 100 e tronco per avere l'indice del bin in cui cade il numero casuale */
      }
      for (int n=0; n<n_bins; n++){
        chi2[i] += pow(V[n]-L/n_bins, 2)/(L/n_bins);  //calcolo il chi-quadro per l'i-esimo blocco
      }
   }

   //Apro il file in cui voglio scrivere i risultati
   std::ofstream file(File_name);

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }

   for(int i=0; i<N; ++i) {
      file << setw(15) << chi2[i] << "\n";
   }
//Chiudo il file
   file.close();
}



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
  andrò a generare un numero M di numeri casuali raggruppati in N blocchi  */
   const int M = 1000000; // Quanti numeri casuali devo generare
   const int N = 100;    // Numero di blocchi in cui li divido
   const int L = M / N;  // Numero di elementi in ogni blocco

   //Genero un numero M di numeri pseudo-casuali
   std::vector<double> r(M);
   for(int i = 0; i < M; ++i) {
      r[i] = rnd.Rannyu();      //random generator uniforme
    }

   //Calcolo la media e la dev std cumulativa
   string File_name = "text1.txt";
   Ave_printer(N, M, L, File_name, 0, r);

   //Calcolo della varianza con errore
   File_name = "text2.txt";
   Ave_printer(N, M, L, File_name, 0.5, r);
   
   //Calcolo dei valori di chi-quadro per gli N blocchi
   int n_bins=100;
   File_name = "text3.txt";
   Chi_quadro(N, L, n_bins, File_name, r);
   
   rnd.SaveSeed();
   return 0;
}