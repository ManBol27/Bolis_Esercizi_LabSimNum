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

int main (int argc, char *argv[]){

  auto start = std::chrono::high_resolution_clock::now(); //voglio stampare il tempo impiegato per runnare la simulazione
  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  double act_temperature = 0.;  //iniziallizzo variabile temperatura attuale
  

  //Apro il file in cui stampare la temperatura al variare dei timesteps
  ofstream temperature_file;
  temperature_file.open("../OUTPUT/temperature_VS_timestep.dat");
  if(!temperature_file.is_open()){
    cerr << "PROBLEM: Unable to open temperature_VS_timestep.dat" << endl;
    exit(EXIT_FAILURE);
  }
  temperature_file << "#   BLOCK:  TEMPERATURE:" << endl;
  double stable_temperature = 0.; //variabile che conterrà la media degli ultimi 100 valori di temperatura del primo blocco
  int ave_steps = 200; //voglio calcolare la temperatura media di equilibrio negli ultimi 200 steps del blocco
  


  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(i == 0){       /*devo farlo per il primo blocco, è quello che 
        mostra un transiente iniziale prima di arrivare ad una temperatura stabile */
        act_temperature = SYS.get_temperature_measurement(); //Get actual temperature (at this timestep)
        temperature_file << setw(12) << i << setw(12) << act_temperature << endl;
        if (j> SYS.get_nsteps()-(ave_steps+1)) stable_temperature += act_temperature; //calcolo la media degli ultimi 200 valori di temperatura
      }
      if(j%10 == 0){
//        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  temperature_file.close();
  stable_temperature /= 200.; //calcolo la media degli ultimi 100 valori di temperatura
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Equilibrium temperature (T*_equilibrium) after "<< SYS.get_nsteps() << " steps computed as the mean temperature of the last " << ave_steps <<" steps is: " << stable_temperature << endl;
  coutf.close();
  SYS.finalize();

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
