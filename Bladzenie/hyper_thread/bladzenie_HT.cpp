#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>
#include <thread>
#include <mutex>
#define NO_threads 8

using namespace std;

mutex m;
fstream plik;
void save_file (int N, double W)
{
    m.lock();

    plik << log10(N) << "   " << log10(W) << endl;

    m.unlock();
}

double variance (double* drunks ,const int K){
    double sum = 0;
    double square_sum = 0;
    for (int i = 0; i < K; i++){
        sum += drunks[i];
        square_sum += (drunks[i]*drunks[i]);
    }
    return (square_sum/(double)K)-pow(sum/(double)K,2);
}

void set_point (int start, int POINTS, const int K , int N){
    for (int i = start; i < POINTS; i++){
    double drunks [K];
    for (int p = 0; p < K; p++){
        int x_coord = 0;
        int random = 0;
        for (int q = 0; q < N; q++){
            random = (rand()%4);
            if (random  < 2)
                x_coord ++;
            else
                x_coord --;
        }
        drunks[p] = (double)x_coord;
    }
    save_file(N , variance(drunks,K));
    N += 10;

    }
}

int main (){

    srand (time(NULL)); // inicjalizacja zarodka

    plik.open("data.txt", ios::out);
    if (plik.good() == false)
        return 0;

    const int K = 10000; // liczba pijanych studentow po zaliczeniu sesji
    int N = 10000; //  pierwotna liczba kroków do akademika
    const int N_POINTS = 10000; // liczba punktów na wykresie

    thread watki[NO_threads];
    int points_start =0;
    int points_stop = N_POINTS/NO_threads;
    int steps_begin = N;

    for (int i=0; i < NO_threads; i++){
        watki[i] = thread (set_point, points_start , points_stop , K , steps_begin);
        points_start = points_stop + 1;
        points_stop += N_POINTS/NO_threads;
        steps_begin = steps_begin + 10*(N_POINTS/NO_threads);
    }


    for (int i=0; i < NO_threads; i++){
        watki[i].join();
    }

    plik.close();

    return 0;
}
