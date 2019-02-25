#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
using namespace  std;

int main (){

    const int max_molecule = 900; // końcowa liczba czastek w symulacji
    const int m_size = 30; // rozmiar układu
    int no_molecules = 1; // początkowa ilość cząstek
    int next[m_size], previous[m_size]; // tablice najbliższych sasiadów
    int matrix[m_size][m_size] = {0}; // nasza macierz
    const int no_probes = 1000; // liczba niezależnych symulacji
    const int MCS = 500; // monte carlo steps


    fstream file;

    file.open("data.txt", ios::out);

    if (file.good() == false)
        cout << "File error" << endl;

    srand(time(NULL));

    for (int i=0; i < m_size; i++){ // uzupełnianie tablicy najbliższych sąsiadów
        next[i] = i+1;
        previous[i] = i -1;
    }
    next[m_size-1] = 0;
    previous[0] = m_size -1;

    double concentration;

    while (no_molecules < max_molecule){

        int x[no_molecules], y[no_molecules]; // tablica położeń cząsteczek
        double diff_coef[no_probes];
        double diff_coefficient;
        double coeff_temp[MCS]; // współczynniki przy każdym kroku czasowym
        concentration = static_cast<double>(no_molecules)/static_cast<double>(m_size*m_size);

        for (int v=0; v < no_probes; v++){

            for (int i=0; i < m_size; i++){ // zerowanie naszej macierzy
                for (int q=0; q < m_size; q++){
                    matrix[i][q] = 0;
                }
            }

            int counter = 0; // licznik cząstek

            while (counter < no_molecules){ // generuję początkowe położenia atomów w macierzy
                int x_temp = (rand()% m_size);
                int y_temp = (rand()% m_size);
                if (matrix[x_temp][y_temp] == 0){
                    x[counter] = x_temp;
                    y[counter] = y_temp;
                    matrix[x_temp][y_temp] = 1;
                    counter ++;
                }
            }

            int random_number =0; // przechowuje wylosowaną liczbę
            int delta_x[no_molecules] = {0};
            int delta_y[no_molecules] = {0};
            int delta_R2[no_molecules];

            for (int i=0; i < MCS; i++){ // kolejne kroki czasowe

                for (int q=0; q < no_molecules; q++){
                    random_number = rand()%4;
                    switch (random_number){
                    case 0:
                        if (matrix[x[q]][next[y[q]]] == 0){
                            matrix[x[q]][y[q]]=0;
                            y[q] = next[y[q]]; // incrementing y coordinate
                            delta_y[q]++;
                            matrix[x[q]][y[q]] = 1; // filling a particle into the matrix
                        }
                        break;
                    case 1:
                        if (matrix[x[q]][previous[y[q]]] == 0){
                            matrix[x[q]][y[q]]=0;
                            y[q] = previous[y[q]]; // decrementing y coordinate
                            delta_y[q]--;
                            matrix[x[q]][y[q]] = 1; // filling a particle into the matrix
                        }
                        break;
                    case 2:
                        if (matrix[next[x[q]]][y[q]] == 0){
                            matrix[x[q]][y[q]]=0;
                            x[q] = next[x[q]]; // incrementing x coordinate
                            delta_x[q]++;
                            matrix[x[q]][y[q]] = 1; // filling a particle into the matrix
                        }
                        break;
                    case 3:
                        if (matrix[previous[x[q]]][y[q]] == 0){
                            matrix[x[q]][y[q]]=0;
                            x[q] = previous[x[q]]; // decrementing x coordinate
                            delta_x[q]--;
                            matrix[x[q]][y[q]] = 1; // filling a particle into the matrix
                        }
                        break;
                   }
                }
                int R2_sum = 0; // sum of R squares, useful for <R2>
                for (int i = 0; i < no_molecules; i++){
                delta_R2[i] = delta_x[i]*delta_x[i] + delta_y[i]*delta_y[i];
                R2_sum += delta_R2[i];
                }
                double R2_avarage = R2_sum/static_cast<double>(no_molecules);
                if (i != 0)
                    coeff_temp[i] = R2_avarage/(4*(i));
                else
                    coeff_temp[i] = 0;
                //cout << coeff_temp[i] << endl;
            }
            double temp_sum =0;
            for (int q=0; q<MCS; q++){
                temp_sum += coeff_temp[q];
            }
            diff_coef[v] = temp_sum/static_cast<double>(MCS);

        } // koniec próby
        double sum_coef = 0; // suma wpó³czynników dla jedego MCS
        for (int i=0; i<no_probes; i++){
            sum_coef +=  diff_coef[i];
        }
        diff_coefficient = sum_coef/(double)no_probes;
        file << concentration << "    " << diff_coefficient << endl;

        no_molecules++;

    }


    file.close();

    return 0;
}
