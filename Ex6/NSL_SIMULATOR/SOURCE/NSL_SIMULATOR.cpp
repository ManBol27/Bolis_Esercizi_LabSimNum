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
#include "system.h"
#include <random>
#include <chrono>

using namespace std;

//Funzione per modificare la temperatura nel file di input
void modifica_temp(double nuova_temp) {
  std::ifstream infile("../INPUT/input.dat");
  std::ofstream outfile("../INPUT/input_mod.dat");  // Scrive in un nuovo file

  if (!infile.is_open() || !outfile.is_open()) {
      std::cerr << "Errore nell'aprire i file!\n";
      return;
  }

  std::string line;
  while (std::getline(infile, line)) {  //legge il file riga per riga e la salva in line
      std::istringstream iss(line);  //istringstream è un oggetto che permette di leggere da una stringa come se fosse un flusso di input (come il "cin")
      std::string keyword;
      iss >> keyword;   //legge la prima parola della riga e la salva in keyword

      // Controlla se la riga contiene la parola chiave "TEMP", allora la modifica con la nuova temperatura
      if (keyword == "TEMP") {
          outfile << "TEMP                   " << nuova_temp << "\n";
      } 
      //altrimenti ristampa la riga così come è
      else {
          outfile << line << "\n";
      }
  }

  infile.close();
  outfile.close();

  // Sovrascrive input.dat con il modificato 
  std::remove("../INPUT/input.dat");
  std::rename("../INPUT/input_mod.dat", "../INPUT/input.dat");
}

int main (int argc, char *argv[]){

  auto start = std::chrono::high_resolution_clock::now(); //voglio stampare il tempo impiegato per runnare la simulazione
  int nconf = 1;
  double current_temperature = 0.5; //temperatura iniziale

  //Iniziallizzo i files in cui verranno stampate le misure al variare della temperatura
  System SYS;
  SYS.initialize_properties_T_(); //Inizializzo i files per le misure al variare della temperatura

  //Creo un ciclo che farà girare NSL_SIMULATOR per ogni temperatura
  for(int i=0; i<16; i++){
    current_temperature = 0.5 + 0.1*i; //sovrascrivo la temperatura caricata, vado da 0.5 a 2 a passi di 0.1
    cout << "Temperatura: " << current_temperature << endl;
    modifica_temp(current_temperature); //modifico il file di input con la nuova temperatura

    //Eseguo il programma NSL_SIMULATOR
    System SYS;
    SYS.initialize();
    SYS.initialize_properties();
    SYS.block_reset(0);
    double act_temperature = 0.;  //iniziallizzo variabile temperatura attuale
  


    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
  //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
     }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    ofstream coutf;
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf.close();
    SYS.finalize();

  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;

  // Stampa il tempo impiegato in secondi
  std::cout << "Tempo impiegato: " << diff.count() << " s" << std::endl;

  return 0;
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
