#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>

using namespace std;

double variance (double* drunks ,const int K){
    double sum = 0;
    double square_sum = 0;
    for (int i = 0; i < K; i++){
        sum += drunks[i];
        square_sum += (drunks[i]*drunks[i]);
    }
    return (square_sum/(double)K)-pow(sum/(double)K,2);
}

int main (){

srand (time(NULL)); // inicjalizacja zarodka

fstream plik;
plik.open("data.txt", ios::out);
if (plik.good() == false)
    return 0;

const int K = 10000; // liczba pijanych studentow po zaliczeniu sesji
int N = 10000; //  pierwotna liczba kroków do akademika
const int N_POINTS = 10000; // liczba punktów na wykresie

for (int i = 0; i < N_POINTS; i++){
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
    plik << log10(N) << "   " << log10(variance(drunks,K)) << endl;
    N += 10;

}

plik.close();

return 0;
}
