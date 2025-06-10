//=========================================================
//                  SIMULATED ANNEALING
//=========================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "../GenRand/random.h"

using namespace std;

// Calcola errore standard cumulativo
double error(const vector<double>& AV, const vector<double>& AV2, int n) {
    if (n == 0) return 0.0;
    return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

// Media cumulativa e errore (ultimo blocco) per blocking
vector<double> Ave_printer(int N, int extractions, int L, const vector<double>& r) {
    vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N);
    for (int i = 0; i < N; ++i) {
        double sum1 = 0.0;
        for (int j = 0; j < L; ++j) sum1 += r[j + i * L];
        ave[i] = sum1 / L;
        av2[i] = ave[i] * ave[i];
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        su2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, su2_prog, i);
    }
    return { sum_prog.back(), err_prog.back() };
}

// Stampa medie cumulative ed errori su file, la uso per stampare energia al variare dei blocchi (fissati i parametri ottimali mu e sigma)
void Ave_printer_blocks(int N, int M, int L, const string& File_name, const vector<double>& r) {
    vector<double> ave(N), av2(N), sum_prog(N), su2_prog(N), err_prog(N);
    for (int i = 0; i < N; ++i) {
        double sum1 = 0.0;
        for (int j = 0; j < L; ++j) sum1 += r[j + i * L];
        ave[i] = sum1 / L;
        av2[i] = ave[i] * ave[i];
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        su2_prog[i] /= (i + 1);
        err_prog[i] = error(sum_prog, su2_prog, i);
    }
    ofstream file(File_name);
    if (!file.is_open()) { cerr << "Error: unable to open " << File_name << endl; exit(1); }
    file << setw(15) << "Media" << setw(15) << "Errore" << "\n";
    for (int i = 0; i < N; ++i) file << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << "\n";
}

// Potenziale V(x)
double V(double x) {
    return pow(x, 4) - 2.5 * pow(x, 2);
}

// Funzione d'onda campione Psi_T
double Psi_T(double x, double mu, double sigma) {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma))
         + exp(-pow(x + mu, 2) / (2 * sigma * sigma));
}

// g(x) per energia locale
double g_x(double x, double mu, double sigma) {
    double t1 = (x - mu) / sigma;
    double t2 = (x + mu) / sigma;
    double num = t1*t1 * exp(-t1*t1/2) + t2*t2 * exp(-t2*t2/2);
    double psi = Psi_T(x, mu, sigma);
    return (1.0/(2*sigma*sigma))*(1 - num/psi) + V(x);
}

// PDF per Metropolis: Psi_T^2
double pdf(double x, double, double, double mu, double sigma) {
    double p = Psi_T(x, mu, sigma);
    return p * p;
}

// Simulated Annealing per ottimizzazione mu,sigma
vector<double> Sim_Annealing(double (*pdf)(double,double,double,double,double), int N, int L, int M, int extractions,
                             const string& process, double mu, double sigma, double T0, Random& rnd) {
    // Apro il file subito, per assicurare che venga creato
    ofstream file("Sim_Anneal_path.dat");
    if (!file.is_open()) { cerr << "Cannot open output file 'Sim_Anneal_path.dat'.\n"; exit(1); }
    cout << "[DEBUG] Opened Sim_Anneal_path.dat for writing\n";
    file << setw(15) << "Mu" << setw(15) << "Sigma" << setw(15) << "Energia"
         << setw(15) << "Err_Energia" << setw(15) << "Temp" << setw(15) << "Passo Metropolis" << "\n";
    file.flush();

    // --- campionamento iniziale ---
    double S = 0.1;
    bool init_ok = false;
    vector<vector<double>> samp;
    for (double s = S; s <= 4.0; s += 0.1) {
        samp = rnd.Metropolis_sampling(pdf, s, M, extractions, process, mu, sigma);
        if (!samp.empty()) { init_ok = true; cout << "Initial sampling OK at S="<<S<<"\n"; S = s; break; }
    }
    if (!init_ok) {
        cerr << "Initial acceptance rate not reachable (S>4).\n";
        file.close();
        return {mu, sigma};
    }

    // energia iniziale
    vector<double> gvals(extractions);
    for (int i = 0; i < extractions; ++i) gvals[i] = g_x(samp[i][0], mu, sigma);
    auto res = Ave_printer(N, extractions, L, gvals);
    double oldE = res[0], oldErr = res[1];
    file << setw(15) << mu << setw(15) << sigma << setw(15) << oldE << setw(15) << oldErr << setw(15) << T0 << setw(15) << S << "\n";
    file.flush(); //serve per stampare direttamente sul file anzichè tenere in memoria nel buffer, per evitare che se ci sono dei bug i risultati già trovati non vengano stampati

    int steps_fixed = 200;
    double S2;

    for (double T = T0; T > 2e-3; T *= 0.97) {
        double delta = max(0.001, 0.02 * (T / T0));
        for (int i = 0; i < steps_fixed; ++i) {
            double mu_p = mu + rnd.Rannyu(-delta, delta);
            double sig_p = sigma + rnd.Rannyu(-delta, delta);
            // campionamento per proposta
            S2 = 0.01;
            bool ok = false;
            vector<vector<double>> s2;
            for (double s = S2; s <= 4.0; s += 0.1) {
                s2 = rnd.Metropolis_sampling(pdf, s, M, extractions, process, mu_p, sig_p);
                if (!s2.empty()) { ok = true; S2 = s; break; }
            }
            if (!ok) continue;  //se ok non è stato impostato a true, allora il codice salta il resto del codice
            for (int k = 0; k < extractions; ++k) {
                gvals[k] = g_x(s2[k][0], mu_p, sig_p);
            }
            res = Ave_printer(N, extractions, L, gvals);
            double E_p = res[0], Err_p = res[1];

            double P = min(1.0, exp(-(E_p - oldE) / T));
            if (rnd.Rannyu() < P) {
                mu = mu_p; sigma = sig_p;
                oldE = E_p; oldErr = Err_p;
            }
        }
        // un solo risultato per temperatura
        file << setw(15) << mu << setw(15) << sigma << setw(15) << oldE
             << setw(15) << oldErr << setw(15) << T << setw(15) << S2 <<"\n";
        file.flush();
    }
    file.close();
    return {mu, sigma};
}

int main() {
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("../GenRand/Primes"); Primes >> p1 >> p2; Primes.close();
    ifstream input("../GenRand/seed.in");
    string prop;
    while (input >> prop)
        if (prop == "RANDOMSEED") input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);

    int M = 5e3, extractions = 1e4;
    const int N = 100;
    int L = extractions / N;
    double T0 = 0.12, mu = 0.5, sigma = 0.7;
    
    
    //Chiamo la funzione per il simulated annealing che restituisce gli ultimi due valori di mu e sigma
    auto opt = Sim_Annealing(pdf, N, L, M, extractions, "Sim_Anneal", mu, sigma, T0, rnd);
    cout << "Last mu, sigma: " << opt[0] << ", " << opt[1] << endl;
    

    // Parte 2: dati per energia media con data blocking (fissati i mu e sigma ottimali)

    //Dal file che ho appena creato leggo le righe ed estraggo i valori mu e sigma che danno minore energia
    std::ifstream results_file ("Sim_Anneal_path.dat");
    if (!results_file.is_open()) {
        std::cerr << "Errore nell'apertura del file.\n";
        return 1;
    }

    double col1, col2, col3, col4, col5, col6;
    double min_energy = 1.e6;
    double opt_mu = 0.;
    double opt_sigma = 0.;
    double opt_S = 0.; 

    string line;

    std::getline(results_file, line); //Leggo la prima linea che contiene le intestazioni
    while (std::getline(results_file, line)) {
        std::istringstream iss(line);
        iss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6;

        if(col3 < min_energy) {
            min_energy = col3;
            opt_mu = col1;
            opt_sigma = col2;
            opt_S = col6;
        }
    }


    cout << "Optimal mu, sigma: " << opt_mu << ", " << opt_sigma << endl;
    cout << "Minimum energy: " << min_energy << endl;
    //Con questi valori di mu e sigma ottimali effettuo l'analisi seguente
    vector<vector<double>> finalSamp;
    finalSamp = rnd.Metropolis_sampling(pdf, opt_S, M, extractions, "Prod", opt_mu, opt_sigma);
        
    vector<double> gvals(extractions);
    for (int i = 0; i < extractions; ++i) gvals[i] = g_x(finalSamp[i][0], opt_mu, opt_sigma);
    Ave_printer_blocks(N, extractions, extractions/N, "Energy_optimal.dat", gvals);


    //Parte 3: campiono dei valori da |psi|^2 da confrontare con soluzioni analitiche
    //NB:qui voglio piu estrazioni per tracciare meglio psi, anche sopra sarebbe stato meglio aumentare questo numero, ma il codice diventa troppo pesante

    int n_samplings = 1e7;
    vector<vector<double>> finalSamp2;
    finalSamp2 = rnd.Metropolis_sampling(pdf, opt_S, M, n_samplings, "Prod", opt_mu, opt_sigma);
    //Stampo su un file i valori estratti da |psi|^2
    std::ofstream psi_file ("Psi_sampled.dat");
   
    for (int i=0; i<n_samplings; i++) {
        psi_file << finalSamp2[i][0] << endl;
    }



    rnd.SaveSeed();
    return 0;
}
