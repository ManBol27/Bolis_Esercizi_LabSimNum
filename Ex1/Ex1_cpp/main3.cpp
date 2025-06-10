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

/* Funzione per la medie cumulative di pigreco con relativo errore */

void Ave_printer(int N, int M, int L, const std::string& File_name, const std::vector<double>& r, double D, double L_ago) {
   /* 
   Iniziallizzo dei vettori che saranno rispettivamente:
     -Vettore dei numeri casuali, vettore delle medie degli N blocchi, idem delle      medie al quadrato, vettore in cui l'elemento i-esimo è la media sui primi i       blocchi, lo stesso per le medie al quadrato e infine un vettore per la dev       std cumulativa 
   */
   std::vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N); 

   
   for(int i = 0; i < N; ++i) {
     double sum1 = 0.0;
     for(int j = 0; j < L; ++j) {
      int k = j + i * L;
      sum1 += r[k];   
   }
     ave[i] = 2*L_ago*L / (sum1*D);    //valore di pigreco restituito dal blocco i-esimo, sum1/L è il numero di hitten versus totale lanci
     av2[i] = pow(ave[i], 2);      //valore di pigreco al quadrato restituito dal blocco i-esimo
       }

   for(int i = 0; i < N; ++i) {
     for(int j = 0; j <= i; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
   }
     sum_prog[i] /= (i + 1);      //media dei valori di pigreco sui primi i blocchi
     su2_prog[i] /= (i + 1);      //media dei valori di pigreco^2 sui primi i blocchi
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

/******************FUNZIONE INTERSECA**********/
/* Funzione che restituisce 0 se il lancio casuale non ha intersecato la linea, 1 se l'ha intersecata*/

double interseca(Random& rnd, double L, double D) {

   double inter = 0;  
   bool incerchio =false;
   double xa = rnd.Rannyu(L, L+D); //primo estremo del segmento (ago lanciato), suppongo WLOG ya=0
   double xb = 0, yb = 0;
   while(!incerchio) {   /*serve per far si che il secondo estremo sia estratto dentro al cerchio inscritto nel quadrato
      vedi jupyter*/
      xb = rnd.Rannyu(xa-L, xa+L); //secondo estremo del segmento (ago lanciato)
      yb = rnd.Rannyu(-L, L); //seconda coordinata del secondo estremo del segmento (ago lanciato)
      if (sqrt(pow(xb-xa,2)+pow(yb,2))<=L) {
         incerchio = true;
      }
   }
   /* Rinormalizzo le coordinate del punto B affinchè l'ago sia lungo L */
   double theta = atan2(yb, xb-xa);
   // double yb_norm = L*sin(theta) ;  non serve per altri scopi
   double xb_norm = xa+L*cos(theta);
   /* Controllo se l'ago interseca la linea */
   if (xb_norm <= L or xb_norm >= L+D) {
      inter = 1;
   }
   return inter;
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

   double L_ago = 10.; //lunghezza ago 
   double D = 30.; //distanza tra le righe 
   //Genero un numero M di estrazioni del tipo intersecato (=1) e non intersecato (=0)

   std::ofstream file1("estrazioni.txt");
   std::vector<double> r_inter(M);
   for(int i = 0; i < M; ++i) {
      r_inter[i] = interseca(rnd, L_ago, D);
      file1 << r_inter[i] << "\n";      
   }

   
   //Calcolo la media di pigreco cumulativa con relativo errore
   string File_name = "simPi.txt";
   Ave_printer(N, M, L, File_name, r_inter, D, L_ago);
   
   rnd.SaveSeed();
   return 0;
}