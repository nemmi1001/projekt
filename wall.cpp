#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "wall.h"


using namespace std;

Wall::Wall(double x_0, double y_0, double x_end, double y_end, double dp){
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));
	

    int n = l/dp + 1;
	cout << "l/dp = " << l/dp << endl;
	double zb_ = l/dp  - (n-1);
	cout << "n = " << n << endl;
	cout << "zb_ = " << zb_ << endl;
	double deltadp = zb_ / (n-1);
	//dp += deltadp;
	cout << "deltadp = " << deltadp << endl;
	
/*	
	if (round(n) != n) {
		printf("\n[Error] Pocet castic neni cele cislo: %d\n", n);
        printf("X0 = [%lf, %lf] a XEND = [%lf, %lf]\n", x_0, y_0, x_end, y_end);
		exit(1);
	}
*/
    this->n = n;

	static int idx_wall = 0;
	idx = idx_wall;
	
	idx_wall++;

	s[0] = (x_end-x_0) / l;
    s[1] = (y_end-y_0) / l;

	double zb_2 = l - sqrt((pow((n-1)*dp*s[0], 2) + pow((n-1)*dp*s[1], 2)));
	cout << "zb_2 = " << zb_2 << endl;
	double deltadp2 = zb_2 / (n-1);
	cout << "deltadp2 =" << deltadp2 << endl;
	dp += deltadp2;
	cout << "dp =" << dp << endl;
	cout << "wall[" << idx << "] upraveno dp = " << dp << endl << endl;
	

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
/*
	double zb = l - (n-1);

	if (zb >= dp/2) {
		double l_ = sqrt((pow(x_end-P[n-1].x, 2) + pow(y_end-P[n-1].y, 2))); 
		P.resize(n+1);
		P[n].x  = P[n-1].x + l_ * s[0];
		P[n].y  = P[n-1].y + l_ * s[1];
		P[n].nx = -s[1];
		P[n].ny =  s[0];
	}
	*/
}

void Wall::create(double x_0, double y_0, double x_end, double y_end, double dp){
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));

    int n = l/dp + 1;
	//cout << "n = " << n;
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
	double zb = l - (n-1);
	if (zb >= dp/2) {
		double l_ = sqrt((pow(x_end-P[n-1].x, 2) + pow(y_end-P[n-1].y, 2))); 
		P.resize(n+1);
		P[n].x  = P[n-1].x + l_ * s[0];
		P[n].y  = P[n-1].x + l_ * s[1];
		P[n].nx = -s[1];
		P[n].ny =  s[0];
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



void solveClosedWallCorners(vector<Wall>& wall) {
	int n = wall.size();

	for (int i = 0; i < n-1; i++)
		solveWallCorner(wall[i], wall[i+1]);
	
	solveWallCorner(wall[n-1], wall[0]);
}

void solveWallCorner(Wall& w1, Wall& w2) {
	double eps = 1e-14;
	
	int n1 = w1.P.size();
	int n2 = w1.P.size();

	if (compareFloatNumbers(w1.P[0].x, w2.P[0].x, eps) && 
		compareFloatNumbers(w1.P[0].y, w2.P[0].y, eps)) {
		
		sumUnitVector(w1.P[0], w2.P[0]);
	}

	else if (compareFloatNumbers(w1.P[n1-1].x, w2.P[n2-1].x, eps) && 
			 compareFloatNumbers(w1.P[n1-1].y, w2.P[n2-1].y, eps)) {

		sumUnitVector(w1.P[n1-1], w2.P[n2-2]);
	}

	else if (compareFloatNumbers(w1.P[0].x, w2.P[n2-1].x, eps) && 
			 compareFloatNumbers(w1.P[0].y, w2.P[n2-1].y, eps)) {
		
		sumUnitVector(w1.P[0], w2.P[n2-2]);
	}

	else if (compareFloatNumbers(w1.P[n1-1].x, w2.P[0].x, eps) && 
			 compareFloatNumbers(w1.P[n1-1].y, w2.P[0].y, eps)) {

		sumUnitVector(w1.P[n1-1], w2.P[0]);
	}

	else {
		cout << "wall[" << w1.idx << "] a wall[" <<  w2.idx << "] nesdílí žádné společné částice" << endl;
	}
}


void defineRectangle(vector<Wall>& wall, double x_0, double y_0, double a, double b, double dp) {
	wall.resize(4);

	wall[0].create(x_0    , y_0	   , x_0    , y_0 + b, dp);
	wall[1].create(x_0    , y_0 + b, x_0 + a, y_0 + b, dp);
	wall[2].create(x_0 + a, y_0 + b, x_0 + a, y_0    , dp);
	wall[3].create(x_0 + a, y_0    , x_0    , y_0	 , dp);

	solveClosedWallCorners(wall);
}

void defineCircle(vector<Wall>& wall, double x_0, double y_0, double r, double dp) {
	vector<double> t;
	double zb;
	linspaceInterval(t, 0, dp/r, 2*M_PI, &zb);

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

	if (zb >= dp/2) {
		
		t.resize(n+1);
		
		wall.resize(n+1);

		wall[n-1].P.resize(2);
		wall[n-1].P[0].x  = x_0 + r*cos(t[n-1]);
		wall[n-1].P[0].y  = y_0 + r*sin(t[n-1]);
		wall[n-1].P[0].nx = cos(t[n-1]);
		wall[n-1].P[0].ny = sin(t[n-1]);

		double l = sqrt((pow(wall[0].P[0].x-wall[n-1].P[0].x, 2) + pow(wall[0].P[0].y-wall[n-1].P[0].y, 2))); 
		double alpha = 2*asin(l/(2*r));
		double s = alpha *r ;
		cout << "s = " << s/4 << endl;
		cout << "dp = " << dp/2 << endl;
		t[n] = t[n-1] + s/(r*2);
		wall[n-1].P[1].x  = x_0 + r*cos(t[n]);
		wall[n-1].P[1].y  = y_0 + r*sin(t[n]);
		wall[n-1].P[1].nx = cos(t[n]);
		wall[n-1].P[1].ny = sin(t[n]);

		wall[n].P.resize(2);
		wall[n].P[0].x  = x_0 + r*cos(t[n]);
		wall[n].P[0].y  = y_0 + r*sin(t[n]);
		wall[n].P[1].x  = x_0 + r*cos(t[0]);
		wall[n].P[1].y  = y_0 + r*sin(t[0]);
		wall[n].P[0].nx = cos(t[n]);
		wall[n].P[0].ny = sin(t[n]);
		wall[n].P[1].nx = cos(t[0]);
		wall[n].P[1].ny = sin(t[0]);

	}

	else {
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
	
}

void defineCircleArc(vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double x_s, double y_s, double dp) {
	vector <double> t;	
	double r      = sqrt((x_0-x_s) * (x_0-x_s) + (y_0-y_s) * (y_0-y_s));
	double l_half = 0.5*sqrt((-x_0+x_end) * (-x_0+x_end) + (-y_0+y_end) * (-y_0+y_end));
	double tau 	  = 2*asin(l_half/r);
	double psi;
	double zb;

	if(x_s == x_0 && y_0 > y_s)
		psi = 0.5*M_PI; 	//horní pi/2

	else if (x_s == x_0 && y_0 < y_s)
		psi = 1.5*M_PI; 	//dolní pi/2

	else {
		double psi_p = atan(fabs(y_0-y_s) / fabs(x_0-x_s));
		
		if(x_s < x_0 && y_s <= y_0)
			psi = psi_p;			//1. kvadrant

		else if(x_s >= x_0 && y_s < y_0)
			psi = M_PI - psi_p; 	//2. kvadrant

		else if(x_s > x_0 && y_s >= y_0)
			psi = M_PI + psi_p;		//3. kvadrant

		else if(x_s <= x_0 && y_s > y_0)
			psi = 2*M_PI - psi_p;	//4. kvadrant
	}
	
	double phi = psi + tau;	
	//cout<<psi<< "   "<< phi <<endl;
	linspaceInterval(t, psi, dp/r, phi, &zb);
	int n = t.size();
	
	
	wall.resize(n);
		for(int i = 0; i < n-1; i++){
		//cout<<t[i]<<endl;
		wall[i].P.resize(2);
		wall[i].P[0].x  = x_s + r*cos(t[i]); 
		wall[i].P[0].y  = y_s + r*sin(t[i]);
		wall[i].P[1].x  = x_s + r*cos(t[i+1]);
		wall[i].P[1].y  = y_s + r*sin(t[i+1]);
		wall[i].P[0].nx = cos(t[i]); 
		wall[i].P[0].ny = sin(t[i]);
		wall[i].P[1].nx = cos(t[i+1]);
		wall[i].P[1].ny = sin(t[i+1]);

	}

    wall[n-1].P.resize(2);
	wall[n-1].P[0].x = wall[n-2].P[1].x;
	wall[n-1].P[0].y = wall[n-2].P[1].y;
	wall[n-1].P[1].x = x_end;
	wall[n-1].P[1].y = y_end;

	wall[n-1].P[0].nx = wall[n-2].P[1].nx;
	wall[n-1].P[0].ny = wall[n-2].P[1].ny;
	wall[n-1].P[1].nx = cos(phi);
	wall[n-1].P[1].ny = sin(phi);	
}


void WallFinalize(vector<Wall>& wall) {
	int n = wall.size();

	for (int i = 0; i < n-1; i++) {

	}
}

void WallSave(vector<Wall>& wall) {
	for (size_t i = 0; i < wall.size(); i++)
		wall[i].save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");
}


void linspaceInterval(vector<double>& t, double t_0, double dt, double t_end, double* zb) {
	int n = (t_end - t_0)/dt;
	*zb = (t_end - t_0)/dt - n;

	cout << "pocet dilku =" <<(t_end - t_0)/dt << endl;
	cout << "zbytek = " << *zb << endl;
	//if (round(n) != n)
	//	cout << "Error[linspace], interval nelze rozedělit "
	t.resize(n);

	for (int i = 0; i < n; i++)
		t[i] = t_0 + i * dt; 

}

bool compareFloatNumbers(double x, double y, double eps) {
	return fabs(x-y) <= eps;
}

void sumUnitVector(particle& P1, particle& P2){
	double nx, ny, l;

	nx = P1.nx + P2.nx;
	ny = P1.ny + P2.ny;
	l  = sqrt(pow(nx, 2) + pow(ny, 2));

	P1.nx = P2.nx = nx / l;
	P1.ny = P2.ny = ny / l;
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