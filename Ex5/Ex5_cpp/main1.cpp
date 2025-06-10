//=========================================================
//                  METROPOLIS CON UNIFORM TRANSITION
//=========================================================

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

/* Funzione per la medie cumulative con dev std cumulativa per il blocking delle estrazioni */

void Ave_printer(int N, int extractions, int L, const std::string& File_name, const std::vector<double>& r) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += r[k];   
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

/*************************************/


//========================================================
//       PDF DA USARE PER METROPOLIS 
//========================================================

//PDF Ground State (è la funzione d'onda al quadrato), scritta in unità del raggio di Bohr (a_0=1)
double p_GroundState(double x, double y, double z) {
  return exp(-2*sqrt(x*x+y*y+z*z))/M_PI;
}

//PDF 2P State (è la funzione d'onda al quadrato), scritta in unità del raggio di Bohr (a_0=1)
double p_2Pstate(double x, double y, double z) {
  return (1./64.)*(2./M_PI)*(x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z))*pow((z/sqrt(x*x+y*y+z*z)),2);
}

//========================================================

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

   // Inizializzo variabili del mio problema: 
   
   /*L'algoritmo di Metropolis implementato nel random.cpp si presenta nella seguente forma 
   Metropolis_sampling (double (*p)(double, double, double), double S, int M, int extractions)
   */

   double S = 1.28; // modulo del passo del random walk
   int M = 1; // numero di step random walk per equilibrare
   int extractions = 1E4; // numero di step dopo l'equilibrazione
   const int N = 100; // numero di blocchi
   int L = extractions / N; // numero di step per ogni blocco
   string process = "Null";  // serve da passare a metropolis per dare un nome specifico al file con acceptance rate 
   //////////////////////////////////////////////////////
   //CALCOLO R MEDIO DI EQUILIBRAZIONE PER GROUND STATE
   cout<<"Procedo con equilibrazione per Ground State"<<endl;
   process = "Equilib_GS";
   std::vector<vector<double>> Equilibration_GS = rnd.Metropolis_sampling (p_GroundState, S, M, extractions, process); // eseguo il random walk
   //Trasformo le posizioni in distanze
   std::vector<double> r_Equilibration_GS (extractions);
   for (int i=0; i<extractions; i++) {
      r_Equilibration_GS[i] = sqrt(pow(Equilibration_GS[i][0],2) + pow(Equilibration_GS[i][1],2) + pow(Equilibration_GS[i][2],2));
   }
   
   std::string File_name = "GS_Equilibration_"+std::to_string(S)+".txt";
   Ave_printer(N, extractions, L, File_name, r_Equilibration_GS);
   cout<<"Processo terminato"<<endl;

   //////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////
  //CALCOLO R MEDIO DOPO EQUILIBRAZIONE PER STATO GS
  cout<<"Procedo con calcolo r_medio Ground State dopo equilibrazione"<<endl;

  M = 1E4;
  extractions = 1E6;
  L = extractions / N;
  process = "GS";
  std::vector<vector<double>> GS = rnd.Metropolis_sampling (p_GroundState, S, M, extractions, process); // eseguo il random walk
   //Trasformo le posizioni in distanze
  std::vector<double> r_GS (extractions);
  for (int i=0; i<extractions; i++) {
      r_GS[i] = sqrt(pow(GS[i][0],2) + pow(GS[i][1],2) + pow(GS[i][2],2));
   }
  File_name = "GS_"+std::to_string(S)+".txt";
  Ave_printer(N, extractions, L, File_name, r_GS);
  cout<<"Processo terminato"<<endl;

   //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //CALCOLO R MEDIO DI EQUILIBRAZIONE PER STATO 2P
  S = 3.; // modulo del passo del random walk
  M = 1; // numero di step random walk per equilibrare
  extractions = 1E4; // numero di step dopo l'equilibrazione
  L = extractions / N; // numero di step per ogni blocco
  cout<<"Procedo con equilibrazione per Stato 2P"<<endl;
  process = "Equilib_2P";
  std::vector<vector<double>> Equilibration_2P = rnd.Metropolis_sampling (p_2Pstate, S, M, extractions, process); 
  //Trasformo le posizioni in distanze
  std::vector<double> r_Equilibration_2P (extractions);
  for (int i=0; i<extractions; i++) {
     r_Equilibration_2P[i] = sqrt(pow(Equilibration_2P[i][0],2) + pow(Equilibration_2P[i][1],2) + pow(Equilibration_2P[i][2],2));
  }
  
  File_name = "2P_Equilibration_"+std::to_string(S)+".txt";
  Ave_printer(N, extractions, L, File_name, r_Equilibration_2P);
  cout<<"Processo terminato"<<endl;

  //////////////////////////////////////////////////////
 
 //////////////////////////////////////////////////////
 //CALCOLO R MEDIO DOPO EQUILIBRAZIONE PER STATO 2P
 cout<<"Procedo con calcolo r_medio stato 2P dopo equilibrazione"<<endl;

 M = 1E4;
 extractions = 1E6;
 L = extractions / N;
 process = "2P";
 std::vector<vector<double>> State_2P = rnd.Metropolis_sampling (p_2Pstate, S, M, extractions, process); // eseguo il random walk
  //Trasformo le posizioni in distanze
 std::vector<double> r_2P (extractions);
 for (int i=0; i<extractions; i++) {
     r_2P[i] = sqrt(pow(State_2P[i][0],2) + pow(State_2P[i][1],2) + pow(State_2P[i][2],2));
  }
 File_name = "2P_"+std::to_string(S)+".txt";
 Ave_printer(N, extractions, L, File_name, r_2P);
 cout<<"Processo terminato"<<endl;

  //////////////////////////////////////////////////////
 //////////////////////////////////////////////////////
  
   
   rnd.SaveSeed();
   return 0;
}