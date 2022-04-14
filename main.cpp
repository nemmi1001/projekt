#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

typedef struct{
	double x, y, nx, ny;
} particle;

void wall(double x_0, double y_0, double x_end, double y_end, double dp); // double ort - ort orientace normály
void save(vector<particle>& P);


int main(){

	wall (2, 2, 5, 6, 1);
	return 0;
}


void wall(double x_0, double y_0, double x_end, double y_end, double dp){ //double ort
	vector<double> x, y, nx, ny;	
	
	double l  = sqrt((pow(x_end-x_0,2) + pow(y_end-y_0,2)));
	double nn = l/dp + 1;
	vector<particle> P(nn);
	double smer[] = {(x_end-x_0) / l, (y_end-y_0) / l};

	P[0].x  = x_0;
	P[0].y  = y_0;
	P[0].nx = -smer[1];
	P[0].ny =  smer[0];

	if (round(nn) != nn) {
		printf("\n[Error] Pocet castic neni cele cislo: %lf\n", nn);
		exit(1);
	}
	else	
		for(int i = 1; i < nn; i++) {
			P[i].x  = P[i-1].x + dp * smer[0];
			P[i].y  = P[i-1].y + dp * smer[1];
			P[i].nx = -smer[1];
			P[i].ny =  smer[0];
		}
	
	save(P);
}

void save(vector<particle>& P) {
	ofstream vystup("wall.txt");

	for(int i = 0; i < P.size(); i++) {
		vystup << fixed << setprecision(2);
		vystup << P[i].x << "\t" << P[i].y << "\t" << P[i].nx << "\t" << P[i].ny << endl;
	}

	vystup.close();
}