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
  

  //Apro il file in cui stampare energia potenziale al variare dei timesteps
  ofstream penergy_file;
  penergy_file.open("../OUTPUT/potential_VS_timestep.dat");
  if(!penergy_file.is_open()){
    cerr << "PROBLEM: Unable to open temperature_VS_timestep.dat" << endl;
    exit(EXIT_FAILURE);
  }
  penergy_file << "#   STEP:  TEMPERATURE:" << endl;
  double act_penergy = 0.; //iniziallizzo variabile energia potenziale attuale


  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(i == 0){       /*devo farlo per il primo blocco, è quello che 
        mostra un transiente iniziale  */
        act_penergy = SYS.get_potential_measurement(); //Get actual penergy (at this timestep)
        penergy_file << setw(12) << i << setw(12) << act_penergy << endl;
      }
      if(j%10 == 0){
//        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  penergy_file.close();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
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
