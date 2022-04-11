#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

void wall(double x_0, double y_0, double x_end, double y_end, double dp); // double ort - ort orientace normály
void save(vector<double>& x, vector<double>& y);



int main(){

	wall (2, 2, 5, 6, 1);
	return 0;
}


void wall(double x_0, double y_0, double x_end, double y_end, double dp){ //double ort
	vector<double> x, y;	

	double l  = sqrt((pow(x_end-x_0,2) + pow(y_end-y_0,2)));
	double nn = l/dp;
	double smer[] = {(x_end-x_0) / l, (y_end-y_0) / l};

	x.push_back(x_0);
	y.push_back(y_0);

	if (round(nn) != nn) {
		printf("\n[Error] Pocet castic neni cele cislo: %lf\n", nn);
		exit(1);
	}
	else	
		for(int i = 0; i < nn; i++) {
			x.push_back(x[i] + dp * smer[0]);
			y.push_back(y[i] + dp * smer[1]);
		}

	save(x, y);
}

void save(vector<double>& x, vector<double>& y) {
	ofstream vystup("wall.txt");

	for(int i = 0; i < x.size(); i++) {
		vystup << fixed << setprecision(2);
		vystup << x[i] << "\t" << y[i] << endl;
	}

	vystup.close();

}