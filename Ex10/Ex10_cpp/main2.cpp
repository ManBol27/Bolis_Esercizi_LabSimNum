//===============================================================================================
// ALGORITMO GENETICO APPLICATO A TRAVELING SALESMAN PROBLEM - PARALLEL COMPUTING (CON MIGRAZONE)
//===============================================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "mpi.h"
#include "../GenRand/random.h"

using namespace std;


//===========================================================================================================================//

// Funzione che scambia una coppia della sequenza, tranne il primo e l'ultimo elemento. G è il numero di geni nell'individuo
vector<int> switcher (int G, vector<int> individuo, Random& rnd) {
   //Scelgo a random due indici del vettore da scambiare tra di loro
   int first_index = static_cast<int>(rnd.Rannyu()*(G-1)) + 1; //prende un indice tra 1 e G-1 (ricorda che l'ultimo indice del vettore è pari a G)
   int sec_index = static_cast<int>(rnd.Rannyu()*(G-1)) + 1; 

   while (sec_index == first_index) {     //evito di scambiare elemento con se stesso
        sec_index = static_cast<int>(rnd.Rannyu() * (G - 1)) + 1;
    }

   swap(individuo[first_index], individuo[sec_index]);
   

   return individuo;
}

//===========================================================================================================================//



//===========================================================================================================================//

// Creo funzione che genera la popolazione iniziale...creo una sequenza ordinata da 0 a G (numero di città o "geni")  che è il primo individuo e 
// gli indiviui successivi sono una permutazione di quello iniziale, a formare una popolazione di M individui
vector<vector<int>> Population_creator (int M, int G, Random& rnd) {

   vector<vector<int>> Population (M, vector<int>(G+1));    //Ogni riga è una sequenza lunga G+1 siccome l'ultima città deve essere la prima

   //Creo tutti gli individui come sequenze ordinate da 0 a G-1 poi le permuterò...
   for (int j=0; j<M; j++) {
      for (int i = 0; i<G; i++) {
         Population[j][i] = i;
      }
      Population[j][G] = 0;
   }
   
   //Genero i restanti individui come il risultato di un numero random di permutazioni di coppie a partire dal primo individuo
   int switches = 0;  // numero di switch per creare un certo individuo, verrà scelto a random
   for (int j=1; j<M; j++) {

      switches = static_cast<int>(rnd.Rannyu()*G) + 1; //sceglie un numero casuale tra 1 e G di scambi
      for (int s=0; s<switches; s++) {
         Population[j] = switcher(G, Population[j], rnd);
      }

   }

   return Population;
}


//==============================================================================================================================//



//===========================================================================================================================//
//Funzione che controlla se un individuo rispetta i vincoli (ogni indice deve comparire una e una sola volta)
//Per verificare ciò, prendo un individuo, lo riordino per valori crescenti, e se è corretto sarà uguale alla sequenza 0,0,1,2,.......,G-1

bool check_func (vector<int> individuo) {


   //Controllo che primo e ultimo elemento siano 0, perchè la sequenza deve avere il formato 0,.......,0 (il cammino deve iniziare e finire nella prima città)
   if (individuo.front() != 0 or individuo.back() != 0) return false;

   // Ordino il vettore in ordine crescente
   std::sort(individuo.begin(), individuo.end());

   //Comparo con la sequenza corretta
   //controllo il primo elemento che deve essere zero
   if (individuo[0] != 0) return false;

   //controllo gli altri elementi se anche solo uno non coincide allora l'individuo controllato non va bene
   for (int i=0; i<individuo.size()-1; i++) {  
      if (individuo[i+1] != i) return false;
   }

   return true;

}

//===========================================================================================================================//



//===========================================================================================================================//
//Funzione che restituisce la fitness di un individuo, calcolata come L^2(x0, x1, ...., x0)
//Come argomento deve prendere un individuo (sequenza considerata), e il vettore delle posizioni 2D delle città

double fitness (vector<int> individuo, vector<vector<double>> cities) {

   double costo = 0. ;
   for (int i = 0; i<individuo.size()-1; i++) {
      costo += (pow((cities[individuo[i]][0] - cities[individuo[i+1]][0]), 2) + pow((cities[individuo[i]][1] - cities[individuo[i+1]][1]), 2));
   }

   return 1. / costo;  //la fitness è l'inversa del costo
}


//===========================================================================================================================//



//===========================================================================================================================//
// Funzione per riordinare gli individui della popolazione da quello con fitness più alta a quello con fitness più bassa

void sort_population_by_fitness(vector<vector<int>>& Population,  const vector<vector<double>>& cities) {
    std::sort(
        Population.begin(), Population.end(), [&](const std::vector<int>& A, const std::vector<int>& B) {
            return fitness(A, cities) > fitness(B, cities);
        }
    );
}


//===========================================================================================================================//




//===========================================================================================================================//
//Funzione che seleziona un individuo (cromosoma) della popolazione con probabilità maggiore se esso ha una buona fitness

vector<int> selector (vector<vector<int>> & Population, int M, Random & rnd) {

   //sceglie casualmente un indice tra 0 e M-1 pesando maggiormente gli indici piccoli, funziona perchè lo applico su una pop gia ordinata per fitness decresc
   int index_chosen = static_cast<int>(M*pow(rnd.Rannyu(), 4.5));   
   return Population[index_chosen];
    
}

//===========================================================================================================================//



//===========================================================================================================================//
//Introduco ora una funzione che può applicare diversi tipi di mutazione ad un individuo ognuno con una certa probabilità

vector<int> Mutator (vector<int> individuo, Random & rnd, int G, double prob_mut) {

   //Introduco un ciclo che potrebbe applicare una di 4 tipi di mutazione 
   double prob = prob_mut;

   //Con una probabilità data da prob_mut eseguo una delle 4 operazioni scelta a caso
   if (rnd.Rannyu() <= prob) {

   int random_mutation = static_cast<int>(rnd.Rannyu()*4); //produce a caso un valore tra 0 e 3

   //Mutazione 1 (switcher)
   if (random_mutation == 0) {
      individuo = switcher(G, individuo, rnd);
   }
 
   //Mutazione 2 (shift di n = int(G/5) posizioni)
   if (random_mutation == 1) {
      int n = static_cast<int>(G/5);
      vector<int> new_individuo (G+1);
      new_individuo = individuo ; //creo una copia dell'individuo originale come vettore di supporto
      for (int i=1; i<G-n; i++) {
         individuo[i] = new_individuo[i+n];
      }
      int l = 0; //contatore di supporto
      for (int j=G-n; j<G; j++) {
         individuo[j] = new_individuo[1+l];
         l++;
      }
   }

   //Mutation 3 (permutazione di blocchi, prende a caso due blocchi di dimensione block_dim e li swappa)
   if (random_mutation == 2) {

      //Scelgo la dim del blocco di elementi contigui che sia massimo 6 ed in ogni caso non eccede la metà della dimensione del vettore
      int block_dim = 6;
      while (block_dim > static_cast<int>(G/2)) block_dim--;

      //Genero due posizioni casuali con vincoli opportuni per non finire fuori range del vettore
      int pos1 = 999999;
      int pos2 = 999999;
      while ((pos1 + block_dim) > G) {
         pos1 = static_cast<int>(rnd.Rannyu()*(G-1)) +1 ; //esclude il primo e l'ultimo elemento che devono restare 0
      }

      while ((pos2 + block_dim) > G or abs(pos1-pos2) < block_dim) {   //deve essere anche che pos1 e pos2 siano distanti piu di block_dim
         pos2 = static_cast<int>(rnd.Rannyu()*(G-1)) +1 ; 
      }

      for (int k = 0; k < block_dim; ++k) {
        std::swap(individuo[pos1 + k], individuo[pos2 + k]);
      }

   }



   //Mutation 4 (lascio la prima e l'ultima città visitata fissate, inverto l'ordine delle altre)
   if (random_mutation == 3) {

      vector<int> new_individuo1 (G+1);
      new_individuo1 = individuo ; //creo una copia dell'individuo originale come vettore di supporto
      for (int i=1; i<G-1; i++) {
         individuo[i] = new_individuo1[G-1-i];
      }
   }

}

   return individuo; 
}
//===========================================================================================================================//




//===========================================================================================================================//
// Funzione che implementa il crossover, prendendo in ingresso due individui, ne restituisce due
vector<vector<int>> Crossover (vector<int> individuo1, vector<int> individuo2, Random & rnd, int G, double prob_cross) {

   //Eseguo la funzione solo con una probabilità prob_cross
   if (rnd.Rannyu() <= prob_cross) {

      //Scelgo un punto casuale dove effettuare il crossover (dal terzo indice al terz'ultimo)
      int cross_index = static_cast<int>((rnd.Rannyu()*(G-3))) +2;

      bool already_in = false;

      vector<int> individuo_temp (G+1);
      individuo_temp = individuo1;  //creo una copia di supporto

      for (int i=cross_index; i<G; i++) {
         for (int j=1; j<G; j++) {
            already_in = false;
            for (int k=0; k<i; k++) {
               if (individuo2[j] == individuo1[k]) {
                  already_in = true;
               }
            }
            if (already_in == false) {
               individuo1[i] = individuo2[j];
               break;
            }
         }
      }


      for (int i=cross_index; i<G; i++) {
         for (int j=1; j<G; j++) {
            already_in = false;
            for (int k=0; k<i; k++) {
               if (individuo_temp[j] == individuo2[k]) {
                  already_in = true;
               }
            }
            if (already_in == false) {
               individuo2[i] = individuo_temp[j];
               break;
            }
         }
      }


   }
   
   return {individuo1, individuo2};

}


//===========================================================================================================================//





//===========================================================================================================================//
//Funzione che crea le posizioni 2D di un numero C di città disposte su una circonferenza di raggio 100 (in coord. cartesiane)

vector<vector<double>> Cities_circumf (int G, Random & rnd) {

   double raggio = 100. ; 
   double theta = 0. ;
   vector<vector<double>> cities (G, vector<double>(2));

   for (int i=0; i<G; i++) {
      theta = rnd.Rannyu(0. , 2*M_PI);
      cities[i][0] = raggio * cos(theta);
      cities[i][1] = raggio * sin(theta);
   }

   return cities;
}

//===========================================================================================================================//



//===========================================================================================================================//
//Funzione che crea le posizioni 2D di un numero C di città disposte dentro un quadrato di lato 100 (in coord. cartesiane)

vector<vector<double>> Cities_square (int G, Random & rnd) {

   double lato = 100. ;
   vector<vector<double>> cities (G, vector<double>(2));

   for (int i=0; i<G; i++) {
      cities[i][0] = rnd.Rannyu(0. , lato);
      cities[i][1] = rnd.Rannyu(0. , lato);
   }

   return cities; 
}

//===========================================================================================================================//


//===========================================================================================================================//
//Funzione per leggere delle coordinate da file (città italiane)

vector<vector<double>> City_Loader (string filename) {

   vector<vector<double>> cities;

   ifstream file(filename);

	if (!file.is_open()) {
		cerr << "Errore nell'apertura del file delle città." << endl;
		return {};
	}

   string line;


   int index = 0;

   while(getline(file, line)) {
      double lat, lon;
      istringstream iss(line);
      iss >> lat >> lon;
      cities.push_back({lat, lon});
      index++;
   }

   return cities;

}



//===========================================================================================================================//



//===========================================================================================================================//
//Funzione per dare la distanza totale del percorso migliore nella popolazione
double best_dist (vector<vector<int>> Popolazione, vector<vector<double>> cities) {
   return 1. / fitness(Popolazione[0], cities);
}
//===========================================================================================================================//




//===========================================================================================================================//
//Funzione per dare la media della distanza totale della prima metà nella popolazione
double mean_dist (vector<vector<int>> Popolazione, vector<vector<double>> cities) {

   double mean_temp = 0. ; 
   for (int i=0; i<static_cast<int>(Popolazione.size()/2); i++) {
      mean_temp = mean_temp + 1. / fitness(Popolazione[i], cities);
   }

   return mean_temp/static_cast<int>(Popolazione.size()/2);
}
//===========================================================================================================================//



//===========================================================================================================================//
// FUNZIONE ALGORITMO GENETICO 

void Genetic (int M, int G, vector<vector<double>> cities, int generations, Random  rnd, double prob_mut, double prob_cross, string type_path, int &rank, vector<int> &Cores, int N_migr) {

   //Preparo i files su cui voglio stampare
   // File 1) dove stampo percorso migliore finale (risultato ottenuto dopo un certo numero di generazioni)
   std::ostringstream filename;
   filename << "Path_"
         << std::fixed << std::setprecision(2)
         << prob_mut << "_"
         << prob_cross << "_"
         << type_path << "_Core"
         << rank
         << ".dat";

   std::ofstream file_path(filename.str());

   if (!file_path.is_open()) {
    std::cerr << "Errore nell'apertura del file: " << filename.str() << std::endl;
    return;
   }

   file_path << "X_coord" << setw(15) << "Y_coord" << endl;


   // File 2) dove stampo lunghezza best e lunghezza media al variare delle generazioni
   std::ostringstream filename2;
   filename2 << "Lenght_"
         << std::fixed << std::setprecision(2)
         << prob_mut << "_"
         << prob_cross << "_"
         << type_path <<"_Core"
         << rank
         << ".dat";

   std::ofstream file_lenght(filename2.str());

   if (!file_lenght.is_open()) {
    std::cerr << "Errore nell'apertura del file: " << filename2.str() << std::endl;
    return;
   }

   file_lenght << "Generation" << setw(15) << "Best distance" << setw(15) << "Mean distance" << endl;


   //Creo un insieme di M individui, ciascuno con G geni. In questo caso sono M sequenze di un numero G di città visitate, con vincoli del problema
   vector<vector<int>> Popolazione = Population_creator (M, G, rnd);
   sort_population_by_fitness (Popolazione, cities);

   double best_path = 0. ; 
   double mean_path = 0. ;
   int temps = 0 ;
   int max_temps = 30;
   
   vector<int> individuo_temp (G+1);
   vector<int> individuo_temp2 (G+1);
   vector<int> individuo_scambiato (G+1); //vettore che riceve l'individuo proveniente da un altro core

   vector<vector<int>> new_pop = Popolazione; //creo una copia di supporto della nuova popolazione


   for (int j=0; j<generations; j++) {

      for (int i=0; i<Popolazione.size(); i+=2) {

         //===============  MIGRAZIONE ============
         //Effettuo lo scambio solo se la generazione è multiplo di N_migr, e solo per il primo individuo i=0
         
         if ((j % N_migr == 0) && (i == 0)) {  
            if (rank == 0) {
               for (int sw=0; sw<3; sw++) {   //esegue tre swap nella sequenza di cores
                  Cores = switcher(Cores.size(), Cores, rnd);
               }
            }
            MPI_Bcast(Cores.data(), Cores.size(), MPI_INT,0, MPI_COMM_WORLD);   //tutti i cores ricevono la stessa sequenza per accordarsi su quali coppie effettuare lo scambio
            MPI_Barrier(MPI_COMM_WORLD);  //in modo tale da assicurarsi che tutti i core siano sincronizzati prima di partire con send/receive

            for (int index=1; index<Cores.size()-1; index+=2) {

               //Core che prima invia poi riceve
               if (rank == Cores[index]) {
                  MPI_Send(Popolazione[i].data(), G+1, MPI_INT, Cores[index+1], 0, MPI_COMM_WORLD);  //.data() è il buffer inviato, ovvero si invia il puntatore al primo elemento del vettore
                  MPI_Recv(individuo_scambiato.data(), G+1, MPI_INT, Cores[index+1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Cores[index] è il rank da cui il messaggio arriva
                  Popolazione[i] = individuo_scambiato;
               }
               //Core che prima riceve poi invia
               else if (rank == Cores[index+1]) {
                  MPI_Recv(individuo_scambiato.data(), G+1, MPI_INT, Cores[index], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Cores[index] è il rank da cui il messaggio arriva
                  MPI_Send(Popolazione[i].data(), G+1, MPI_INT, Cores[index], 0, MPI_COMM_WORLD);  //Cores[index] è il rank a cui inviare il messaggio
                  Popolazione[i] = individuo_scambiato;

               }
            }
         }
         //=============== END-MIGRAZIONE =========
         individuo_temp = selector(Popolazione, M, rnd);
         individuo_temp2 = selector(Popolazione, M, rnd);
         temps = 0;
         //Nell'ottica di fare crossover impongo di prendere due genitori diversi
         while (temps < max_temps && individuo_temp == individuo_temp2) {
            individuo_temp2 = selector(Popolazione, M, rnd);
            temps++;
         }
         //Inserisco possibilità di crossover
         auto results = Crossover(individuo_temp, individuo_temp2, rnd, G, prob_cross);
         individuo_temp = results[0];
         individuo_temp2 = results[1];
         //Inserisco possibilità di mutation
         individuo_temp = Mutator(individuo_temp, rnd, G, prob_mut);
         individuo_temp2 = Mutator(individuo_temp2, rnd, G, prob_mut);
         new_pop[i] = individuo_temp;
         new_pop[i+1] = individuo_temp2;
      }

      //Riordino la nuova popolazione in base alla fitness decrescente
      sort_population_by_fitness (new_pop, cities);

      //Calcolo lunghezza del percorso migliore
      best_path = best_dist (new_pop, cities);

      //Calcolo lunghezza media del percorso per la metà della popolazione (individui migliori)
      mean_path = mean_dist (new_pop, cities);

      file_lenght << j+1<< setw(15) << best_path << setw(15) << mean_path << endl;

      Popolazione = new_pop;

   }

   for (int k=0; k<G+1; k++) {
   file_path << cities[Popolazione[0][k]][0] << setw(15) << cities[Popolazione[0][k]][1] << endl;
   }
   
   file_lenght.close();
   file_path.close();

   return;

}


//===========================================================================================================================//



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




   //Prova per vedere se la funzione check va correttamente

   int M = 400; // numero di individui nella popolazione
   int G = 34; //numero di geni per individuo (aka numero di città da visitare)
   int generations = 700 ; //numero di generazioni, ovvero numero di volte che viene creata una nuova popolazione a partire da quella precedente
   int N_migr = 30;
   

   /*
   GENERA PUNTI DISTRIBUITI SU CIRCONFERENZA E DENTRO QUADRATO - EX9
   //Genero le posizioni 2D di G città su una circonferenza con raggio 100
   vector<vector<double>> cities_circ = Cities_circumf(G, rnd);

   //Genero le posizioni 2D di G città entro un quadrato di lato 100
   vector<vector<double>> cities_squar = Cities_square(G, rnd);
   */


   //Carica le coordinate delle città da file
   string file_name = "cap_prov_ita.dat";

   vector<vector<double>> ita_cities = City_Loader(file_name);

   G = ita_cities.size();
   


   //ALGORITMO GENETICO 
   double prob_mut = 0.45;
   double prob_cross = 0.8;


   //INTRODUCO UNA PARALLELIZZAZIONE DEL CODICE
   int size, rank;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);     //la variabile size assume il valore del numero di processi attivati dall'utente al lancio del programma
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);     //la variabile rank assume il numero del core che sta eseguendo questo codice
   if(size>9){cout<<"Hai scelto troppi processi"<<endl;
      return 1;
   }

   if(size%2 != 0){cout<<"Scegli un numero pari di processi"<<endl;
      return 1;
   }


   //GENERO VETTORE DEI CORES del tipo {0,1,2,....size-1}
   //mi serve per generare diverse configurazioni di scambio, permutando l'ordine e prendendo rank adiacenti nel vettore cores, 
   //saranno quelli che si scambiano individui tra di loro
   vector<int> Cores (size);
   for (int i=0; i<size; i++) {
      Cores[i] = i;
   }

   // usa seed diverso per ogni rank, così che partano da popolazioni iniziali differenti
   seed[0] = seed[0] + rank*1;
   seed[1] = seed[1] + rank*2;
   seed[2] = seed[2] + rank*3;
   seed[3] = seed[3] + rank*4;

   rnd.SetRandom(seed, p1, p2);


   /*
   //Per quadrato
   string type_path = "Square";
   Genetic (M, G, cities_squar, generations, rnd, prob_mut, prob_cross, type_path, rank);


   //Resetto i seed del random number generator per avere stabilità dei risultati
	rnd.SetRandom(seed,p1,p2);
	

   //Per cerchio
   type_path = "Circle";
   Genetic (M, G, cities_circ, generations, rnd, prob_mut, prob_cross, type_path, rank);
   */

   //Per città italiane
   string type_path = "CittaIta";
   Genetic (M, G, ita_cities, generations, rnd, prob_mut, prob_cross, type_path, rank, Cores, N_migr);

   MPI_Finalize();
  
   
   rnd.SaveSeed();
   return 0;
}