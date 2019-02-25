#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <mutex>
#include <vector>
#define NO_threads 8

using namespace std;

const double T_max = 3.5; //  maks. temperatura zredukowana
const double T_min = 1.5; // min. temperatura zredukowana
const double T_step = 0.01; // krok temperaturowy
const unsigned int  m_size = 100;
const unsigned int MCS = 230000;
const unsigned int rej_time = 30000; // ilosc krokow na poczatku ktore zostaja odrzucone
const unsigned int probe_time = 100; // co ktora probka jest rozpatrywana
const double T_critical = 2.26929; // temperatura krytyczna

mutex m;

fstream time_log; // informacje o czaise dzialania programu
fstream file_scalling; // skalowanie
fstream file_mag_avarage; // srednia magnetyzacja
fstream file_susceptibility; // podatnosc magnetyczna
fstream file_energy; // energia
fstream file_thermal_conductivity; // pojemnosc cieplna
fstream file_Binder_cumulant; // kumulanty Bindera

void save_to_file (double time, double magnetization, double susceptibility, double energy, double thermal_cond, double Binder, double scal_x, double scal_y){
    m.lock();

    file_mag_avarage << time << "  " << magnetization << endl;
    file_susceptibility << time << "  " << susceptibility << endl;
    file_energy << time << "  " << energy << endl;
    file_thermal_conductivity << time << "  " << thermal_cond << endl;
    file_Binder_cumulant << time << "  " << Binder << endl;
    file_scalling << scal_x << "  " << scal_y << endl;

    m.unlock();
}

double avarage (vector<double> tab){
    double sum = 0;
    for (unsigned int i=0; i < tab.size(); i++){
        sum += tab[i];
    }
    return sum/static_cast<double>(tab.size());
}

double variance (vector<double> tab){
    double sum_square = 0;
    for (unsigned int i = 0; i < tab.size(); i++){
        sum_square += tab[i]*tab[i];
    }
    return (sum_square/static_cast<double>(tab.size())) - (avarage(tab)*avarage(tab));
}

double get_avarage (int table[][m_size]){
    double temp_sum = 0;
    for (unsigned int i = 0; i<m_size; i++){
        for (unsigned int j = 0; j<m_size; j++){
            temp_sum += table[i][j];
        }
    }
    return temp_sum/(m_size*m_size);
}

double get_energy (int table[][m_size], int previous [], int next[]){
    double H_temp = 0;
    for (unsigned int i=0; i < m_size; i++){
        for (unsigned int j=0; j < m_size; j++){
            H_temp -= table[i][j]*(table[next[i]][j] + table[previous[i]][j] + table[i][next[j]] + table[i][previous[j]]);
        }
    }
    return H_temp/static_cast<double>(2);
}

double get_cumulant (vector<double> magnet){
    vector<double> magnet_2;
    vector<double> magnet_4;
    for (unsigned int i=0; i < magnet.size(); i++){
        magnet_2.push_back(magnet[i]*magnet[i]);
        magnet_4.push_back(pow(magnet[i],4));
    }
    return 1-(avarage(magnet_4)/(3*avarage(magnet_2)*avarage(magnet_2)));
}

void set_point(double T_begin, double T_end){

    int matrix[m_size][m_size];
    int next[m_size], previous[m_size];
    const int dU_possible[5] = {-8,-4,0,4,8}; // mo¿liwe wartoœci dU

    for (unsigned int i=0; i < m_size; i++){ // wypelnienie tablicy najblizszych sasiadow
        next[i] = i+1;
        previous[i] = i -1;
    }
    next[m_size-1] = 0;
    previous[0] = m_size -1;

    srand(time(NULL));

    while (T_begin < T_end){ //petla po temperaturach

        vector<double> magnetization;
        vector<double> energy;
        double susceptibility;
        double thermal_conductivity;

        double exp_table[5]; // tworzê tablice mo¿liwych  eksponent
        for (int i=0; i<5; i++){ // uzupe³nianie tablicy mo¿liwych eksponent
            exp_table[i] = exp(-dU_possible[i]/T_begin);
        }

        for (unsigned int i = 0; i < m_size; i++){ //generuje macierz
            for (unsigned int q = 0; q < m_size; q++){
                unsigned int x = rand()%2;
                if (x == 0)
                    matrix[i][q] = -1;
                else
                    matrix[i][q] = 1;
            }
        }

        double dU; //roznica energii
        unsigned int dU_iterator = 0;

        for (unsigned int step = 0; step < MCS; step++){
            for (unsigned int i = 0; i < m_size; i++){
                for (unsigned int j = 0; j < m_size; j++){
                    dU = 2* matrix[i][j]*(matrix[next[i]][j] + matrix[previous[i]][j] + matrix[i][next[j]] + matrix[i][previous[j]]);
                    for (int q = 0; q < 5; q++){
                        if (dU == dU_possible[q])
                            dU_iterator = q;
                    }
                    if (static_cast<double>(rand())/static_cast<double>(RAND_MAX) < min(exp_table[dU_iterator], static_cast<double>(1))){
                        matrix[i][j] = (-matrix[i][j]); // zamieniam spiny
                    }
                }
            }

            if (step > rej_time){
                if (step%probe_time == 0){

                    magnetization.push_back(abs(get_avarage(matrix))); // obliczam magnetyzacje
                    energy.push_back(get_energy(matrix, previous, next)); // obliczam energie
            }
            }
        }

        susceptibility = (variance(magnetization)*m_size*m_size/T_begin);
        thermal_conductivity = variance(energy)/((T_begin*T_begin)*(m_size*m_size));
        double x_scalling = log(abs(1- T_begin/T_critical)) + log(m_size);
        double y_scalling = log(avarage(magnetization)) + 0.125*log(m_size);
        save_to_file(T_begin,avarage(magnetization), susceptibility , avarage(energy)/static_cast<double>(m_size*m_size),
            thermal_conductivity, get_cumulant(magnetization), x_scalling, y_scalling); // zapis do pliku
        T_begin += T_step; // inkrementacja temperatury
    }
}

int main (){

    double time_start = time (NULL);
    file_scalling.open("data_scalling.txt", ios::out);
    file_mag_avarage.open("data_mag_av.txt", ios::out);
    file_susceptibility.open("data_susceptibility.txt", ios::out);
    file_energy.open("data_energy.txt", ios::out);
    file_thermal_conductivity.open("data_thermal.txt", ios::out);
    file_Binder_cumulant.open("data_Binder.txt", ios::out);
    time_log.open("time_log.txt", ios::out);

    thread core[NO_threads];

    double points_start = T_min;
    double points_stop = T_min + (T_max - T_min)/NO_threads;

    for (int i=0; i < NO_threads; i++){ //tworze watki
        core[i] = thread(set_point, points_start, points_stop);
        points_start = points_stop + T_step;
        points_stop += (T_max - T_min)/NO_threads;
    }

    for (int i=0; i < NO_threads; i++){ //czekam na watki
        core[i].join();
    }

    file_scalling.close();
    file_mag_avarage.close();
    file_susceptibility.close();
    file_energy.close();
    file_thermal_conductivity.close();
    file_Binder_cumulant.close();

    time_log << "Czas dzialania programu: " << static_cast<double>(time (NULL)) - time_start << " s" << endl;

    time_log.close();

    return 0;
}
