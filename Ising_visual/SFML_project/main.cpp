#include <SFML/Graphics.hpp>
#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

//const double T = 1.7;// temperatura zredukowana
const unsigned int  m_size = 200;
const unsigned int MCS = 23000000;
const unsigned int rej_time = 3000; // ilosc krokow na poczatku ktore zostaja odrzucone
const unsigned int probe_time = 100; // co ktora probka jest rozpatrywana

void fill_neighbours(int next_el[m_size], int previous_el[m_size]) {
	for (unsigned int i = 0; i < m_size; i++) { // wypelnienie tablicy najblizszych sasiadow
		next_el[i] = i + 1;
		previous_el[i] = i - 1;
	}
	next_el[m_size - 1] = 0;
	previous_el[0] = m_size - 1;

}

void generate_matrix(int MATRIX[m_size][m_size]) {
	for (unsigned int i = 0; i < m_size; i++) { //generuje macierz
		for (unsigned int q = 0; q < m_size; q++) {
			unsigned int x = rand() % 2;
			if (x == 0)
				MATRIX[i][q] = -1;
			else
				MATRIX[i][q] = 1;
		}
	}

}

void draw_matrix(int MATRIX[m_size][m_size], sf::RenderWindow &okno) {
	sf::RectangleShape rectangle(sf::Vector2f(120, 50)); // tworze sobie prostokat do rysowania
	rectangle.setPosition(0, 0);
	for (unsigned int i = 0; i < m_size; i++) {
		for (unsigned int j = 0; j < m_size; j++) {
			if (MATRIX[i][j] == 1) {
				rectangle.setFillColor(sf::Color::Black);
				okno.draw(rectangle);
			}
			else if (MATRIX[i][j] == -1) {
				rectangle.setFillColor(sf::Color::White);
				okno.draw(rectangle);
			}
			if (j != m_size -1)
				rectangle.move(0, 4);
		}
		if (i != m_size -1)
			rectangle.setPosition(4 * (i+1), 0);
	}

	okno.display();
}

int main() {
	
	double T; // temperatura zredukowana
	cout << "Podaj temperature zredukowana: " << endl;
	cin >> T;
	sf::RenderWindow oknoAplikacji(sf::VideoMode(m_size*4, m_size*4, 32), "Ising model");
	oknoAplikacji.clear(sf::Color(0, 0, 0));
	int matrix[m_size][m_size];
	int next[m_size], previous[m_size]; // tablica najblizszych s¹siadów
	const int dU_possible[5] = { -8,-4,0,4,8 }; // mo?liwe warto?ci dU

	fill_neighbours(next, previous); // wype³niam tablice najblizszych s¹siadów

	srand(time(NULL));

	double exp_table[5]; // tworz? tablice mo?liwych  eksponent
	for (int i = 0; i<5; i++) { // uzupe?nianie tablicy mo?liwych eksponent
		exp_table[i] = exp(-dU_possible[i] / T);
	}

	generate_matrix(matrix);

	double dU; //roznica energii
	unsigned int dU_iterator = 0;
	unsigned int step = 0;

	while (oknoAplikacji.isOpen() && (step < MCS))
	{
		
		sf::Event zdarzenie;
		while (oknoAplikacji.pollEvent(zdarzenie))
		{
			if (zdarzenie.type == sf::Event::Closed)
			{
				oknoAplikacji.close();
			}
		}
		for (unsigned int i = 0; i < m_size; i++) {
			for (unsigned int j = 0; j < m_size; j++) {
				dU = 2 * matrix[i][j] * (matrix[next[i]][j] + matrix[previous[i]][j] + matrix[i][next[j]] + matrix[i][previous[j]]);
				for (int q = 0; q < 5; q++) {
					if (dU == dU_possible[q])
						dU_iterator = q;
				}
				if (static_cast<double>(rand()) / static_cast<double>(RAND_MAX) < min(exp_table[dU_iterator], static_cast<double>(1))) {
					matrix[i][j] = (-matrix[i][j]); // zamieniam spiny
				}
			}
		}

		if (step > rej_time) {
			if (step%probe_time == 0) {
				
				draw_matrix(matrix, oknoAplikacji);

			}
		}
		
		step++;
	}
	
	return 0;
}