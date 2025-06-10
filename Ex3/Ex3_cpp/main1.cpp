#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> 
#include "../GenRand/random.h"

using namespace std;

double error(const std::vector<double>& AV, const std::vector<double>& AV2, int n) {
    if (n == 0)
	return 0.0;
    else
	return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

/****************************************************/
/******************** FUNZIONE PUT_CALL_diretto **********/

/* Funzione per la medie cumulative del prezzo delle opzioni PUT e CALL con relativo errore */

void PUT_CALL_diretto(int N, int M, int L, const std::string& File_name, const std::vector<double>& r, double S0, double T_delivery, double K, double r_free, double volatility, double put_call) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle  medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += exp(-r_free*T_delivery)*std::max(0., put_call*((S0*exp((r_free-0.5*pow(volatility, 2))*T_delivery + volatility*r[k]))-K));   // prezzo dell'opzione al tempo attuale
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

   /************ PARAMETRI FINANZIARI DELLA SIMULAZIONE ******************/
   const double S0 = 100.; //prezzo dell'asset sottostante al tempo attuale t0
   const double T_delivery = 1.; // delivery time -> tempo in anni dopo il quale scade l'opzione
   const double K = 100.; //strike price
   const double r_free = 0.1; //tasso di interesse risk-free annualizzato
   const double volatility = 0.25; //volatilità del sottostante (sigma)
   /*********************************************************************/



   //Genero un numero M di numeri pseudo-casuali estratti da una gaussiana N(0, t), in questo caso scelgo t=1 
   std::vector<double> r_gauss01(M);
   for(int i = 0; i < M; ++i) {
      r_gauss01[i] = rnd.Gauss(0,1);      //random generator da distribuzione gaussiana con media=0 e varianza=1
    }

   //Calcolo la media cumulativa del prezzo CALL con relativo errore
   double put_call = 1.; //1 per CALL, -1 per PUT
   string File_name = "CALL_diretto.txt";
   PUT_CALL_diretto(N, M, L, File_name, r_gauss01, S0, T_delivery, K, r_free, volatility, put_call);

   //Calcolo la media cumulativa del prezzo PUT con relativo errore
   put_call = -1.;   //1 per CALL, -1 per PUT
   File_name = "PUT_diretto.txt";
   PUT_CALL_diretto(N, M, L, File_name, r_gauss01, S0, T_delivery, K, r_free, volatility, put_call);


   
   rnd.SaveSeed();
   return 0;
}