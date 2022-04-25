#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "wall.h"


using namespace std;

Wall::Wall(double x_0, double y_0, double x_end, double y_end, double dp){
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));
	
    double n = l/dp + 1;
	if (round(n) != n) {
		printf("\n[Error] Pocet castic neni cele cislo: %lf\n", n);
        printf("X0 = [%lf, %lf] a XEND = [%lf, %lf]\n", x_0, y_0, x_end, y_end);
		exit(1);
	}
    this->n = n;

	static int idx_wall = 0;
	idx = idx_wall;
	idx_wall++;

	s[0] = (x_end-x_0) / l;
    s[1] = (y_end-y_0) / l;

	P.resize(n);
	P[0].x  = x_0;
	P[0].y  = y_0;
	P[0].nx = -s[1];
	P[0].ny =  s[0];
		
    for(int i = 1; i < n; i++) {
        P[i].x  = P[i-1].x + dp * s[0];
        P[i].y  = P[i-1].y + dp * s[1];
        P[i].nx = -s[1];
        P[i].ny =  s[0];
    }
}

void Wall::create(double x_0, double y_0, double x_end, double y_end, double dp){
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));
	//double alpha = 2 * asin(l / (2*r));

	cout << "l = " << l << endl;
    int n = l/dp + 1;
	//cout << "n = " << n << endl; 
	if (round(n) != n) {
		printf("\n[Error] Pocet castic neni cele cislo: %d\n", n);
        printf("X0 = [%lf, %lf] a XEND = [%lf, %lf]\n", x_0, y_0, x_end, y_end);
		exit(1);
	}
    this->n = n;

	static int idx_wall = 0;
	idx = idx_wall;
	idx_wall++;

	s[0] = (x_end-x_0) / l;
    s[1] = (y_end-y_0) / l;

	P.resize(n);
	P[0].x  = x_0;
	P[0].y  = y_0;
	P[0].nx = -s[1];
	P[0].ny =  s[0];
		
    for(int i = 1; i < n; i++) {
        P[i].x  = P[i-1].x + dp * s[0];
        P[i].y  = P[i-1].y + dp * s[1];
        P[i].nx = -s[1];
        P[i].ny =  s[0];
    }
}

void Wall::save(vector<particle>& P, const string& filename) {
	ofstream vystup(filename);

	for(size_t i = 0; i < P.size(); i++) {
		vystup << fixed << setprecision(2);
		vystup << P[i].x << "\t" << P[i].y << "\t" << P[i].nx << "\t" << P[i].ny << endl;
		
	}

	vystup.close();
}



void solveWallCorners(vector<Wall>& wall) {
	WallCorner(wall[0].P, wall[1].P);
	WallCorner(wall[1].P, wall[2].P);
	WallCorner(wall[2].P, wall[3].P);
	WallCorner(wall[3].P, wall[0].P);
}

void WallCorner(vector<particle>& P1, vector<particle>& P2) {
	double nx, ny, l;
	double eps = 1e-14;

	int n1 = P1.size();
	int n2 = P2.size();

	if (compareFloatNumbers(P1[0].x, P2[0].x, eps) && 
		compareFloatNumbers(P1[0].y, P2[0].y, eps)) {

		nx = P1[0].nx + P2[0].nx;
		ny = P1[0].ny + P2[0].ny;
		l  = sqrt(pow(nx, 2)+pow(ny, 2));

		P1[0].nx = P2[0].nx = nx / l;
		P1[0].ny = P2[0].ny = ny / l;
	}

	else if (compareFloatNumbers(P1[n1-1].x, P2[n2-1].x, eps) && 
			 compareFloatNumbers(P1[n1-1].y, P2[n2-1].y, eps)) {

		nx = P1[n1-1].nx + P2[n2-1].nx;
		ny = P1[n1-1].ny + P2[n2-1].ny;
		l  = sqrt(pow(nx, 2)+pow(ny, 2));

		P1[n1-1].nx = P2[n2-1].nx = nx / l;
		P1[n1-1].ny = P2[n2-1].ny = ny / l;
	}

	else if (compareFloatNumbers(P1[0].x, P2[n2-1].x, eps) && 
			 compareFloatNumbers(P1[0].y, P2[n2-1].y, eps)) {

		nx = P1[0].nx + P2[n2-1].nx;
		ny = P1[0].ny + P2[n2-1].ny;
		l  = sqrt(pow(nx, 2)+pow(ny, 2));

		P1[0].nx = P2[n2-1].nx = nx / l;
		P1[0].ny = P2[n2-1].ny = ny / l;
	}

	else if (compareFloatNumbers(P1[n1-1].x, P2[0].x, eps) && 
			 compareFloatNumbers(P1[n1-1].y, P2[0].y, eps)) {
				 
		nx = P1[n1-1].nx + P2[0].nx;
		ny = P1[n1-1].ny + P2[0].ny;
		l  = sqrt(pow(nx, 2)+pow(ny, 2));

		P1[n1-1].nx = P2[0].nx = nx / l;
		P1[n1-1].ny = P2[0].ny = ny / l;
	}

	else {
		cout << "Stěny nesdílí žádné společné body" << endl;
	}
}


void defineRectangle(vector<Wall>& wall, double x_0, double y_0, double a, double b, double dp) {
	wall.resize(4);

	wall[0].create(x_0    , y_0	   , x_0    , y_0 + b, dp);
	wall[1].create(x_0    , y_0 + b, x_0 + a, y_0 + b, dp);
	wall[2].create(x_0 + a, y_0 + b, x_0 + a, y_0    , dp);
	wall[3].create(x_0 + a, y_0    , x_0    , y_0	 , dp);

	solveWallCorners(wall);
}

void defineCircle(vector<Wall>& wall, double x_0, double y_0, double r, double dp) {
	vector<double> t;
	linspace(t, 0, dp/r, 2*M_PI);
	
	int n = t.size();
	wall.resize(n);

	for(int i = 0; i < n-1; i++){
		wall[i].P.resize(2);
		wall[i].P[0].x  = x_0 + r*cos(t[i]); 
		wall[i].P[0].y  = y_0 + r*sin(t[i]);
		wall[i].P[1].x  = x_0 + r*cos(t[i+1]);
		wall[i].P[1].y  = y_0 + r*sin(t[i+1]);
		wall[i].P[0].nx = cos(t[i]); 
		wall[i].P[0].ny = sin(t[i]);
		wall[i].P[1].nx = cos(t[i+1]);
		wall[i].P[1].ny = sin(t[i+1]);

	}

	wall[n-1].P.resize(2);
	wall[n-1].P[0].x  = x_0 + r*cos(t[n-1]);
	wall[n-1].P[0].y  = y_0 + r*sin(t[n-1]);
	wall[n-1].P[1].x  = x_0 + r*cos(t[0]);
	wall[n-1].P[1].y  = y_0 + r*sin(t[0]);
	wall[n-1].P[0].nx = cos(t[n-1]);
	wall[n-1].P[0].ny = sin(t[n-1]);
	wall[n-1].P[1].nx = cos(t[0]);
	wall[n-1].P[1].ny = sin(t[0]); 
}

void defineCirlceArc(vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double r, double dp) {

}


void WallFinalize(vector<Wall>& wall) {

}

void WallSave(vector<Wall>& wall) {
	for (size_t i = 0; i < wall.size(); i++)
		wall[i].save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");
}


void linspace(vector<double>& t, double t_0, double dt, double t_end) {
	int n = (t_end - t_0)/dt; 
	//if (round(n) != n)
	//	cout << "Error[linspace], interval nelze rozedělit "
	t.resize(n);
	t[0] = t_0;

	for (int i = 0; i < n; i++)
		t[i] = i * dt; 

}

bool compareFloatNumbers(double x, double y, double eps) {
	return fabs(x-y) <= eps;
}

double integrateRectMethod(double x, double dx, double f(double x)) {
	double fx = f(x);
	return fx * dx;
}

double integrateTrapMethod(double x, double dx, double f(double x), double a, double b) {
	double fx = f(x);
	if (x==a|| x == b){
		return 1/2. * fx * dx;
	}
	return fx * dx;
}

double integrateFunction(double f(double x), double a, double b){
	int N = 1e5;
	double x, dx = (b-a)/N;
	double S = 0;

	for(int i = 1; i <= N; i++) {
		x = a + i * dx;
		//S += integrateRectMethod(x, dx, f);
		S += integrateTrapMethod(x, dx, f, a, b);
	}
	
	return S;
}