#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#define NO_threads 8

using namespace std;

const double PI = 3.14159265358979323846;

const double E_max = 4.; //  maks. pole elektrostatyczne
const double E_min = 0.; // min. pole eletrostatyczne
const double E_step = 0.01; // krok pola
const unsigned int x_size = 40;
const unsigned int y_size = 20;
const unsigned int MCS = 230000;
const unsigned int rej_time = 30000; // ilosc krokow na poczatku ktore zostaja odrzucone
const unsigned int probe_time = 100; // co ktora probka jest rozpatrywana
const int trial_angle = 10; // probny kat

const double N_e = 1.7;
const double N_0 = 1.5;
const double ksi = 20;

fstream time_log; // informacje o czaise dzialania programu
fstream file_refractive;
fstream matrix_info;

void draw_matrix(int table[x_size][y_size]){
    for (int i=0; i < y_size; i++){
        for (int q=0; q < x_size; q++){
            matrix_info << table[i][q] << " ";
        }
        matrix_info << endl;
    }
    matrix_info << endl;
    cout << "Kuniec" << endl;
}

double rad(int angle){
    return (static_cast<double>(angle)/static_cast<double>(180))*PI;
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

double get_ref_index (int table[x_size][y_size]){
    double ref_sum = 0;
    for (unsigned int i = 0; i<x_size; i++){
        for (unsigned int j = 0; j<y_size; j++){
            ref_sum += N_e*N_0/sqrt(N_0*N_0*pow(cos(rad(table[i][j])),2)+ N_e*N_e*pow(sin(rad(table[i][j])),2));
        }
    }
    return ref_sum/(x_size*y_size);
}


int main (){

    long double time_start = time (NULL);
    file_refractive.open("data_refractive.txt", ios::out);
    time_log.open("time_log.txt", ios::out);

    int matrix[x_size][y_size];
    int next_x[x_size], previous_x[x_size], next_y[x_size], previous_y[y_size]; // tablica najbl. sasiadow

    for (unsigned int i=0; i < x_size; i++){ // wypelnienie tablicy najblizszych sasiadow w X
        next_x[i] = i+1;
        previous_x[i] = i -1;
    }
    next_x[x_size-1] = 0;
    previous_x[0] = x_size -1;

    for (unsigned int i=0; i < y_size; i++){ // wypelnienie tablicy najblizszych sasiadow w Y
        next_y[i] = i+1;
        previous_y[i] = i -1;
    }
    next_y[y_size-1] = 0;
    previous_y[0] = y_size -1;



    double P_table[181]; // tablica wielomianow LG dla cosin
    for (int x = 0; x <= 180; x++){
        P_table[x] = (3./2.)*(cos(rad(x))*cos(rad(x))) -0.5;
    }

    srand(time(NULL));

    double E_begin = E_min;

    while (E_begin < E_max){ //petla po energiach

        vector<double> refractive_index;

        for (unsigned int i = 0; i < x_size; i++){ //generuje macierz
            for (unsigned int q = 0; q < y_size; q++){
                if (q == 0 || q == (y_size-1)){
                    matrix[i][q] = 0;
                }
                else {
                    int x = rand()% 181;
                    if (x <= 90)
                        matrix[i][q] = x;
                    else
                        matrix[i][q] = 90-x;
                }
            }
        }

        double U_old; // energia starej konfiguracji
        double U_new; // energia nowej konfiguracji
        double dU; //roznica energii
        int new_angle; // nowy probny kat

        for (unsigned int step = 0; step < MCS; step++){
            for (unsigned int i = 0; i < x_size; i++){
                for (unsigned int j = 1; j < (y_size-1); j++){
                    U_old = -ksi*(P_table[abs(matrix[i][j]-matrix[i][next_y[j]])]+P_table[abs(matrix[i][j]-matrix[i][previous_y[j]])]
                            +P_table[abs(matrix[i][j]-matrix[next_x[i]][j])]+P_table[abs(matrix[i][j]-matrix[previous_x[i]][j])])-
                            E_begin*E_begin*P_table[abs(90-matrix[i][j])];

                    double rand_number = (static_cast<double>(rand())/RAND_MAX) - 0.5;
                    new_angle = matrix[i][j] + static_cast<int>(rand_number*trial_angle);
                    if (new_angle > 90)
                        new_angle -= 180;
                    else if (new_angle < (-90))
                        new_angle += 180;

                    U_new = -ksi*(P_table[abs(new_angle-matrix[i][next_y[j]])]+P_table[abs(new_angle-matrix[i][previous_y[j]])]
                            +P_table[abs(new_angle-matrix[next_x[i]][j])]+P_table[abs(new_angle-matrix[previous_x[i]][j])])-
                            E_begin*E_begin*P_table[abs(90-new_angle)];

                    dU = U_new - U_old;
                    if (dU < 0)
                        matrix[i][j] = new_angle;
                    else if ((static_cast<double>(rand())/RAND_MAX) < min(exp(-dU), static_cast<double>(1))){
                        matrix[i][j] = new_angle; // zamieniam k¹ty
                    }
                }
            }

            if (step > rej_time){
                if (step%probe_time == 0){

                    refractive_index.push_back(get_ref_index(matrix)); // obliczam sredni wspolczynnik
            }
            }
        }

        file_refractive << E_begin << "  " << avarage(refractive_index) << endl;
        E_begin += E_step; // inkrementacja pola elektrostatycznego
    }


    file_refractive.close();

    time_log << "Czas dzialania programu: " << static_cast<long double>(time (NULL)) - time_start << " s" << endl;

    time_log.close();

    return 0;
}
