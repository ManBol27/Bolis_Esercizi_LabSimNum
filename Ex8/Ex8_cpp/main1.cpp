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

void Ave_printer(int N, int extractions, int L, const std::string& File_name, const std::vector<double>& r, double mu, double sigma) {
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
   std::ofstream file(File_name, std::ios::app);  //apro in append mode

// Controllo la corretta apertura del file
   if (!file.is_open()) {
       cerr << "Error: unable to open " + File_name << endl;
       exit(1);
     }

   //Stampo soltanto il valore complessivo della media cumulativa e della dev std cumulativa
   file << setw(15) << mu << setw(15) << sigma << setw(15) << sum_prog.back() << setw(15) << err_prog.back() << "\n";
   
//Chiudo il file
   file.close();
   }

/**************************************/

/*************************************/


//=========================================================
//       Funzione g(x) di cui calcolare la media


double V(double x) {
    return pow(x, 4) - 2.5 * pow(x, 2);
}

double Psi_T(double x, double mu, double sigma) {
    double term1 = exp(-pow(x - mu, 2) / (2 * pow(sigma, 2)));
    double term2 = exp(-pow(x + mu, 2) / (2 * pow(sigma, 2)));
    return term1 + term2;
}

double g_x(double x, double mu, double sigma) {
    double numerator =
        pow((x - mu) / sigma, 2) * exp(-pow(x - mu, 2) / (2 * pow(sigma, 2))) +
        pow((x + mu) / sigma, 2) * exp(-pow(x + mu, 2) / (2 * pow(sigma, 2)));

    double psi = Psi_T(x, mu, sigma);

    return (1.0 / (2 * pow(sigma,2))) * (1 - (numerator / psi)) + V(x);
}



//=========================================================


//========================================================
//       PDF DA USARE PER METROPOLIS 

double pdf(double x, double y, double z, double mu, double sigma) {
   return Psi_T(x, mu, sigma) * Psi_T( x,  mu, sigma);
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
   
   /*L'algoritmo di Metropolis implementato nel random.cpp AGGIORNATO PER QUESTO ESERCIZIO si presenta nella seguente forma 
   Metropolis_sampling (double (*p)(double, double, double), double S, int M, int extractions, double mu, double sigma)
   */

   double S = 0.1; // modulo del passo del random walk
   int M = 1E4; // numero di step random walk per equilibrare
   int extractions = 1E4; // numero di step dopo l'equilibrazione
   const int N = 100; // numero di blocchi
   int L = extractions / N; // numero di step per ogni blocco
   string process = "Null";  // serve da passare a metropolis per dare un nome specifico al file con acceptance rate 
   
  
  //////////////////////////////////////////////////////
  //CAMPIONAMENTI DA PDF DOPO EQUILIBRAZIONE
  cout<<"Procedo con calcolo r_medio Ground State dopo equilibrazione"<<endl;

  M = 1E4;
  extractions = 1E6;
  L = extractions / N;
  process = "Ex8";

  //Preparo il file che stamperà i valori di g(x) per ogni coppia di parametri (mu, sigma) con relativi errori
  string File_name = "g(x)_mu_sigma.dat";
  //Apro il file in cui voglio scrivere i risultati
  std::ofstream file(File_name);  

  // Controllo la corretta apertura del file
     if (!file.is_open()) {
         cerr << "Error: unable to open " + File_name << endl;
         exit(1);
       }
  
     //Stampo soltanto il valore complessivo della media cumulativa e della dev std cumulativa
     file << setw(15) << "mu" << setw(15) << "sigma" << setw(15) << "Media" << setw(15) << "Errore" << "\n";
     
  //Chiudo il file
     file.close();

    
 // Inizializzo il vettore che conterrà le estrazioni dopo la fase di equilibrazione   
  std::vector<double> g_x_extracted (extractions);
  
  //CICLO SUI MU E SIGMA PER TROVARE I VALORI DI g(x)
  //Per diverse coppie di parametri (mu, sigma) estraggo valori x_i dalla densità di probabilità
  for (double mu=0.0; mu<=1.2; mu+=0.05) {
      for (double sigma=0.0; sigma<=1.2; sigma+=0.05) {
         S = 0.1; // modulo del passo del random walk reiniziallizzato
         std::vector<vector<double>> Campionamento = rnd.Metropolis_sampling (pdf, S, M, extractions, process, mu, sigma); 
         while (Campionamento.empty()) {
            //cout << "Acceptance rate not in range [0.48, 0.52]. Retrying..." << endl;
            Campionamento = rnd.Metropolis_sampling (pdf, S, M, extractions, process, mu, sigma);
            if (S > 4.) {
               cout<<"Per mu = "<<mu<<" e sigma = "<<sigma<<" acceptance rate del 50 '/, non raggiunto"<<endl;
               break;
            }
            S = S+0.1; // Incremento S per cercare di ottenere un acceptance rate corretto
         }
         
         if (S<=4. && !Campionamento.empty()) {

         cout<<"Raggiunto acceptance rate di circa 50'/, per mu = "<<mu<<" e sigma = "<<sigma<<endl;

         // Con i valori estratti da Campionamento, calcolo i corrispettivi valori di g(x)
         g_x_extracted.resize(extractions);
         for (int i=0; i<extractions; i++) {
            g_x_extracted[i] = g_x(Campionamento[i][0], mu, sigma);
         }
         //Facendo la media cumulativa sui g(x) estratti, calcolo il valore di aspettazione dell'hamiltoniana
         // e la sua dev std cumulativa
         Ave_printer(N, extractions, L, File_name, g_x_extracted, mu, sigma);
         }
      }
   }   
   

 cout<<"Processo terminato"<<endl;

  //////////////////////////////////////////////////////
 //////////////////////////////////////////////////////
  
   
   rnd.SaveSeed();
   return 0;
}