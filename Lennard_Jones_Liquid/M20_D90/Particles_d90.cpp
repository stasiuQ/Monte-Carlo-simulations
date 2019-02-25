#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#define NO_threads 8

using namespace std;

const double PI = 3.14159265358979323846;
const double R_interaction = 2.5;
const double x_step = 0.01;
const double y_step = 0.01;
const double r_step = 0.01;
const double dr = 0.01;
const double density = 0.90; // gestosc czastek
const double T = 0.7; // temperatura zredukowana ukladu
const unsigned int m_size = 20;
const unsigned int MCS = 230000;
const unsigned int rej_time = 30000;
const unsigned int probe_time = 100;
const unsigned int no_molecules = static_cast<int>(density*m_size*m_size);
const double Rmax = static_cast<double>(m_size/2.-1.0);

fstream time_log; // informacje o czaise dzialania programu
fstream file_probability;
fstream matrix_info;

void initialize_matrix (double x_pos[no_molecules], double y_pos[no_molecules]){
    const double free_space_x = 0.0;
    const double free_space_y = 0.0;
    double x_counter =0.+ free_space_x;
    double y_counter = 0.+ free_space_y;
    for (unsigned int i=0; i < no_molecules; i++){
        if (x_counter < m_size -free_space_x){
            x_pos[i] = x_counter;
            y_pos[i] = y_counter;
            x_counter += 1.0 + free_space_x;
        }
        else{
            x_counter = 0. + free_space_x;
            y_counter += 1.0 + free_space_y;
            x_pos[i] = x_counter;
            y_pos[i] = y_counter;
        }
    }
}

double Probability(double rad, double* x, double* y){
    double delta_R;
    double dx,dy;
    double temp_probability =0;
    unsigned int counter;
    if (rad < dr)
        return 0;
    for (unsigned int i=0; i < no_molecules; i++){
        counter =0;
        for (unsigned int j = 0; j < no_molecules; j++){
            dx = abs( x[i] - x[j]);
            dy = abs( y[i] - y[j]);
            if(dx > m_size/2.)
                dx = m_size - dx;
            if(dy > m_size/2.)
                dy = m_size - dy;

            delta_R = sqrt(pow(dx,2) + pow(dy,2));
            if ((delta_R >= rad-dr)&&(delta_R <= rad))
                counter ++;
        }
        temp_probability += counter/(PI*2*dr*rad*density);
    }
    return temp_probability/(static_cast<double>(no_molecules));
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

void save_to_file (vector<double>* tab, double* rad_values, unsigned int rad_number){
    for (unsigned int i=0; i < rad_number; i++){
        file_probability << rad_values[i] << "  " << avarage(tab[i]) << endl;
    }
}

void draw_matrix (double* x_tab, double* y_tab){
    for (unsigned int i = 0; i < no_molecules; i++){
        matrix_info << x_tab[i] << "   " << y_tab[i] << endl;
    }
    matrix_info << endl << endl;
}

int main (){

    long double time_start = time (NULL);
    file_probability.open("data_probability.txt", ios::out);
    time_log.open("time_log.txt", ios::out);
    matrix_info.open("matrix_info.txt", ios::out);

    double x_position[no_molecules];
    double y_position[no_molecules];

    initialize_matrix(x_position, y_position);

    srand(time(NULL));

    /*-------TABLICE WEKTOROW PRAWDOPODOBIENSTW------------------*/
    unsigned int radius_number = static_cast<int>(Rmax/r_step);
    double radius_values[radius_number];
    vector<double> probability_values [radius_number];

    double r_counter = 0; // pomocniczy do wypelnienia tablicy promieni
    for (unsigned int i=0; i < radius_number; i++){
        radius_values[i] = r_counter;
        r_counter += r_step;
    }

    /////////////////////////////////////////

    double dU; //roznica energii
    double x_new_position, y_new_position; // probne polozenia
    double dx , dy , R_square; // stare wzgledne polozenia
    double dx_new, dy_new, R_new_square; // nowe wzgledne  polozenia

    /* -----ALGORYTM METROPOLISA------------*/

    for (unsigned int step = 0; step < MCS; step++){
        for (unsigned int i = 0; i < no_molecules; i++){
            dU =0;
            x_new_position = x_position[i] + (static_cast<double>(rand())/RAND_MAX-0.5)*x_step;
            y_new_position = y_position[i] + (static_cast<double>(rand())/RAND_MAX-0.5)*y_step;

            /* sprawdzam periodyczne warunki brzegowe*/

            if(x_new_position > m_size)
                x_new_position -= m_size;
            if(x_new_position < 0. )
                x_new_position += m_size;
            if(y_new_position > m_size)
                y_new_position -= m_size;
            if(y_new_position < 0.)
                y_new_position += m_size;

            for (unsigned int j=0; j < no_molecules; j++){ // pętla po oddzialywujacych czastkach
                if (j != i){
                    dx = abs(x_position[i]-x_position[j]);
                    dy = abs(y_position[i]-y_position[j]);

                    if(dx > m_size/2.)  // sprawdzam warunku brzegowe dla odleglosci
                        dx = m_size - dx;
                    if(dy > m_size/2.)
                        dy = m_size - dy;
                    R_square = dx*dx + dy*dy;

                    dx_new = abs(x_new_position-x_position[j]);
                    dy_new = abs(y_new_position-y_position[j]);

                    if(dx_new > m_size/2.)
                        dx_new = m_size - dx_new;
                    if(dy_new > m_size/2.)
                        dy_new = m_size - dy_new;
                    R_new_square = dx_new*dx_new + dy_new*dy_new;

                    if (R_new_square <= 1.0) // przerywamy pętle kiedy czastki by sie pokrywaly
                        break;

                    if ((R_square <= R_interaction) && (R_new_square <= R_interaction)){ // jezeli czastki oddzialywuja
                        dU += 4.*(1./pow(R_new_square,6) - 1./pow(R_square,6));
                    }
                }
            }
            if((R_new_square > 1.0) && (static_cast<double>(rand())/RAND_MAX <= min(1.0,exp(-dU/T)))){ // pierwszy wyraz jest poprawka po przerwaniu (break)
                x_position[i] = x_new_position;
                y_position[i] = y_new_position;
            }
        }
        /*-----KONIEC PETLI PO CZASTKACH-------*/

        if (step > rej_time){
            if (step%probe_time == 0){

                for (unsigned int i=0; i < radius_number; i++){
                    probability_values[i].push_back(Probability(radius_values[i], x_position, y_position));
                }
                cout << (static_cast<double>(step)/MCS)*100. << " %" << endl;
        }
        }
    }
    /*------------KONIEC ALG. METROPOLISA------------*/


    save_to_file(probability_values, radius_values, radius_number);
    draw_matrix(x_position, y_position);


    file_probability.close();
    matrix_info.close();

    time_log << "Czas dzialania programu: " << static_cast<long double>(time (NULL)) - time_start << " s" << endl;

    time_log.close();

    return 0;
}
