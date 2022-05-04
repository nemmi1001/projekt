#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "wall.h"


using namespace std;

Wall::Wall(double x_0, double y_0, double x_end, double y_end, double dp) {
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));

    int n = l/dp + 1;
	this->n = n;

	s[0] = (x_end-x_0) / l;
    s[1] = (y_end-y_0) / l;

	P.resize(n);
	P[0].x  = x_0;
	P[0].y  = y_0;
	P[0].nx = -s[1];
	P[0].ny =  s[0];
		
    for (int i = 1; i < n; i++) {
        P[i].x  = P[i-1].x + dp * s[0];
        P[i].y  = P[i-1].y + dp * s[1];
        P[i].nx = -s[1];
        P[i].ny =  s[0];
    }

	double zb = sqrt((pow(x_end-P[n-1].x, 2) + pow(y_end-P[n-1].y, 2)));

	if (zb > 1e-14) {
		cout << endl << "Note: wall nelze rozdelit rovnomerne" << endl;
		cout << "posledni particle posunuta o dp = " << zb << endl << endl;

		P.resize(n+1);
		this->n++;

		P[n].x  = P[n-1].x + zb * s[0];
		P[n].y  = P[n-1].y + zb * s[1];
		P[n].nx = -s[1];
		P[n].ny =  s[0];
	}
}

void Wall::create(double x_0, double y_0, double x_end, double y_end, double dp) {
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));

    int n = l/dp + 1;
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
		
    for (int i = 1; i < n; i++) {
        P[i].x  = P[i-1].x + dp * s[0];
        P[i].y  = P[i-1].y + dp * s[1];
        P[i].nx = -s[1];
        P[i].ny =  s[0];
    }

	double zb = sqrt((pow(x_end-P[n-1].x, 2) + pow(y_end-P[n-1].y, 2)));

	if (zb > 1e-14) {
		cout << endl << "Note: wall[" << idx << "] nelze rozdelit rovnomerne s dp = " << dp << endl;
		cout << "posledni particle posunuta o dp = " << zb << endl << endl;

		P.resize(n+1);
		this->n++;

		P[n].x  = P[n-1].x + zb * s[0];
		P[n].y  = P[n-1].y + zb * s[1];
		P[n].nx = -s[1];
		P[n].ny =  s[0];
	}
}

void Wall::create_fit(double x_0, double y_0, double x_end, double y_end, double dp) {
	l = sqrt((pow(x_end-x_0, 2) + pow(y_end-y_0, 2)));

    int n = l/dp + 1;
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

	double zb = l - sqrt((pow((n-1)*dp*s[0], 2) + pow((n-1)*dp*s[1], 2)));
	if (zb > 1e-14) {
		double dpold = dp;
		double deltadp = zb / (n-1);
		dp += deltadp;
		cout << endl << "Note: wall[" << idx << "] nelze rozdelit rovnomerne s dp = " << dpold << endl;
		cout << "wall[" << idx << "] rozdelena rovnomerne s dp = " << dp << endl << endl;
	}
		
    for (int i = 1; i < n; i++) {
        P[i].x  = P[i-1].x + dp * s[0];
        P[i].y  = P[i-1].y + dp * s[1];
        P[i].nx = -s[1];
        P[i].ny =  s[0];
    }
}

void Wall::save(vector<particle>& P, const string filename) {
	ofstream vystup(filename);

	for (size_t i = 0; i < P.size(); i++) {
		vystup << fixed << setprecision(4);
		vystup << P[i].x << "\t" << P[i].y << "\t" << P[i].nx << "\t" << P[i].ny << endl;
		
	}

	vystup.close();
}

void Wall::save2VTK(vector<particle>& P, const string filename) {
int np = P.size();
        
ofstream vystup(filename);

// file header
vystup << "# vtk DataFile Version 2.0" << endl;
vystup << "SPH boundary points and normals" << endl;
vystup << "ASCII" << endl;
vystup << "DATASET POLYDATA" << endl;

// souradnice bodu a jejich normaly
vystup << "POINTS " << np << " float" << endl;
for (int i = 0; i < np; i++)
	vystup << P[i].x << " " << P[i].y << " " << 0 << endl;

vystup << endl << "POINT_DATA " << np << endl;
vystup << "NORMALS " << "normaly" << " float" << endl;
for (int i = 0; i < np; i++)
	vystup << P[i].nx << " " << P[i].ny << " " << 0 << endl;

/* --------------------- BARVY ------------------------- */
vystup << endl << "SCALARS barva_bodu float 1" << endl;
vystup << "LOOKUP_TABLE tabulka_barev" << endl;
for (int i = 0; i < np; i++)
	vystup << "0.0" << endl;
	
vystup << endl << "SCALARS barva_normal float 1" << endl;
vystup << "LOOKUP_TABLE tabulka_barev" << endl;
vystup << "0.0" << endl;
for (int i = 0; i < np-1; i++)
	vystup << "1.0" << endl;

vystup << endl << "LOOKUP_TABLE tabulka_barev " << 2 << endl;
vystup << "0.0 0.5 0.0 1.0" << endl; // zelená
vystup << "1.0 0.0 0.0 1.0 "<< endl; // červená

vystup.close();
}


void solveClosedWallCorners(vector<Wall>& wall) {
	int n = wall.size();

	for (int i = 0; i < n-1; i++)
		solveWallCorner(wall[i], wall[i+1]);
	
	solveWallCorner(wall[n-1], wall[0]);
}

void solveNWallCorners(vector<Wall>& wall) {
	int n = wall.size();

	for (int i = 0; i < n-1; i++)
		solveWallCorner(wall[i], wall[i+1]);
}

void solveWallCorner(Wall& w1, Wall& w2) {
	double eps = 1e-14;
	
	int n1 = w1.P.size();
	int n2 = w1.P.size();

	if (compareFloatNumbers(w1.P[0].x, w2.P[0].x, eps) && 
		compareFloatNumbers(w1.P[0].y, w2.P[0].y, eps)) {

		cout << "Spatna orientace sten wall[" << w1.idx << "] a wall[" << w2.idx << "]" << endl;
		sumUnitVector(w1.P[0], w2.P[0]);
	}

	else if (compareFloatNumbers(w1.P[n1-1].x, w2.P[n2-1].x, eps) && 
			 compareFloatNumbers(w1.P[n1-1].y, w2.P[n2-1].y, eps)) {
				 
		cout << "Spatna orientace sten wall[" << w1.idx << "] a wall[" << w2.idx << "]" << endl;
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

	else
		cout << endl << "Note: wall[" << w1.idx << "] a wall[" <<  w2.idx << "] nesdílí žádné společné částice" << endl << endl;
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
	// rozdeli interval presne podle parametru dp, vlozi posledni particle podle toho jestli 
	// je vzdalena alespon o dp/2r od prvni castice
	
	vector<double> t;
	double zb;
	
	linspaceInterval(t, 0, dp/r, 2*M_PI, zb);
	cout << endl << "Note: wall nelze rozdelit rovnomerne s dp = " << dp << endl;
	cout << "posledni particle posunuta o dp = " << zb * r << endl << endl;

	int n = t.size();
	wall.resize(n);

	for (int i = 0; i < n-1; i++){
		wall[i].idx = i;
		wall[i].n   = 2;
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
	
	wall[n-1].idx = n-1;
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

void defineCircle_fit(vector<Wall>& wall, double x_0, double y_0, double r, double dp) {
	// upravi parametr dp a rozdeli interval rovnomerne

	vector<double> t;

	linspaceInterval_fit(t, 0, dp/r, 2*M_PI, r);
	
	int n = t.size() - 1;
	wall.resize(n);

	for (int i = 0; i < n-1; i++){
		wall[i].idx = i;
		wall[i].n   = 2;
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
	
	wall[n-1].idx = n-1;
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

void defineCircleArc(vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double x_s, double y_s, int ipar, double dp) {
	vector <double> t;	
	double r      = sqrt((x_0-x_s) * (x_0-x_s) + (y_0-y_s) * (y_0-y_s));
	double l_half = 0.5 * sqrt((-x_0+x_end) * (-x_0+x_end) + (-y_0+y_end) * (-y_0+y_end));
	double tau 	  = 2 * asin(l_half/r);
	double psi;
	double zb;

	if (x_s == x_0 && y_0 > y_s)
		psi = 0.5 * M_PI; 	//horní pi/2

	else if (x_s == x_0 && y_0 < y_s)
		psi = 1.5 * M_PI; 	//dolní pi/2

	else {
		double psi_p = atan(fabs(y_0-y_s) / fabs(x_0-x_s));
		
		if(x_s < x_0 && y_s <= y_0)
			psi = psi_p;			//1. kvadrant

		else if(x_s >= x_0 && y_s < y_0)
			psi = M_PI - psi_p; 	//2. kvadrant

		else if(x_s > x_0 && y_s >= y_0)
			psi = M_PI + psi_p;		//3. kvadrant

		else if(x_s <= x_0 && y_s > y_0)
			psi = 2 * M_PI - psi_p;	//4. kvadrant
	}
	
	double phi = psi + tau;	
	linspaceInterval(t, psi, dp/r, phi, zb);
	cout << endl << "Note: wall nelze rozdelit rovnomerne s dp = " << dp << endl;
	cout << "posledni particle posunuta o dp = " << zb * r << endl << endl;
	int n = t.size();
	
	wall.resize(n);

	for (int i = 0; i < n-1; i++){
		wall[i].P.resize(2);
		wall[i].P[0].x  = x_s + r*cos(t[i]); 
		wall[i].P[0].y  = y_s + r*sin(t[i]);
		wall[i].P[1].x  = x_s + r*cos(t[i+1]);
		wall[i].P[1].y  = y_s + r*sin(t[i+1]);
		wall[i].P[0].nx = ipar * cos(t[i]); 
		wall[i].P[0].ny = ipar * sin(t[i]);
		wall[i].P[1].nx = ipar * cos(t[i+1]);
		wall[i].P[1].ny = ipar * sin(t[i+1]);
	}

    wall[n-1].P.resize(2);
	wall[n-1].P[0].x = wall[n-2].P[1].x;
	wall[n-1].P[0].y = wall[n-2].P[1].y;
	wall[n-1].P[1].x = x_end;
	wall[n-1].P[1].y = y_end;

	wall[n-1].P[0].nx = ipar * wall[n-2].P[1].nx;
	wall[n-1].P[0].ny = ipar * wall[n-2].P[1].ny;
	wall[n-1].P[1].nx = ipar * cos(phi);
	wall[n-1].P[1].ny = ipar * sin(phi);	
}

void defineCircleArcAlt(vector<Wall>& wall, double x_0, double y_0, double x_b, double y_b, double r, int ori, int ipar, double dp) {
double x_a    = 0.5 * (x_b+x_0);	// průsečík osy úsečky s úsečkou
double y_a    = 0.5 * (y_0+y_b);  
double l_half = sqrt((x_a-x_0) * (x_a-x_0) + (y_a-y_0) * (y_a-y_0));
double h      = sqrt(r*r-l_half * l_half);
			
			
double x_1 = (y_a * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
 			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			 - y_0 * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			 + h*h*x_0 - h*h*x_a - r*r*x_0 + r*r*x_a - x_0*x_a*x_a - x_0*x_0*x_a + x_0*y_0*y_0 
			 + x_0*y_a*y_a + x_a*y_0*y_0 + x_a*y_a*y_a + x_0*x_0*x_0 + x_a*x_a*x_a - 2*x_0*y_0*y_a 
			 - 2*x_a*y_0*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double x_2 = (y_0 * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a))
			 - y_a * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a))
			 + h*h*x_0 - h*h*x_a - r*r*x_0 + r*r*x_a - x_0*x_a*x_a - x_0*x_0*x_a + x_0*y_0*y_0 + x_0*y_a*y_a 
			 + x_a*y_0*y_0 + x_a*y_a*y_a + x_0*x_0*x_0 + x_a*x_a*x_a - 2*x_0*y_0*y_a - 2*x_a*y_0*y_a)
			 / (2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double y_1 = (x_0 * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			 - x_a * sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			 * (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			 + h*h*y_0 - h*h*y_a - r*r*y_0 + r*r*y_a + x_0*x_0*y_0 + x_0*x_0*y_a + x_a*x_a*y_0 
			 + x_a*x_a*y_a - y_0*y_a*y_a - y_0*y_0*y_a + y_0*y_0*y_0 + y_a*y_a*y_a - 2*x_0*x_a*y_0 
			 - 2*x_0*x_a*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double y_2 = (x_a*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			* (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			- x_0*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)
			* (-h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) 
			+ h*h*y_0 - h*h*y_a - r*r*y_0 + r*r*y_a + x_0*x_0*y_0 + x_0*x_0*y_a + x_a*x_a*y_0 + x_a*x_a*y_a 
			- y_0*y_a*y_a - y_0*y_0*y_a + y_0*y_0*y_0 + y_a*y_a*y_a - 2*x_0*x_a*y_0 - 2*x_0*x_a*y_a)
			/ (2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));
		
//cout << x_1 << "\t" << y_1 << endl;
//cout << x_2 << "\t" << y_2 << endl;

// ori= 1  konvex, 
// ori=-1  konkáv

double px_konvex, py_konvex;
double px_konkav, py_konkav;
double p_x, p_y;

px_konvex = py_konvex = px_konkav = py_konkav = p_x = p_y = 0;

if (y_1 > y_2) {
	py_konvex = y_1;
	px_konvex = x_1;

	py_konkav = y_2;
	px_konkav = x_2;
}

else if (y_1 <= y_2) {
	py_konvex = y_2;
	px_konvex = x_2;

	py_konkav = y_1;
	px_konkav = x_1;
}


if (ori == 1) {
	p_x = px_konvex;
	p_y = py_konvex;
}

else if (ori == -1){
	p_x = px_konkav;
	p_y = py_konkav;
}

else
	cout << "Špatná orientace zadej 1 nebo -1" << endl;

defineCircleArc(wall, x_0, y_0, x_b , y_b,  p_x,  p_y, ipar, dp);
}


void saveEachWall(vector<Wall>& wall) {
	for (size_t i = 0; i < wall.size(); i++)
		wall[i].save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");
}

void saveEachBoundary(vector<Wall>& Bndr) {
	for (size_t i = 0; i < Bndr.size(); i++)
		Bndr[i].save(Bndr[i].P, "bndr" + std::to_string(i + 1) + ".txt");
}


void createPartBoundary(Wall& B, vector<Wall>& wall) {
    int i, j, np, nbp, nw = wall.size();

	static int idx_bndr = 0;
	B.idx = idx_bndr;
	idx_bndr++;

    for (i = 0; i < nw; i++) {
        np  = wall[i].P.size();
        nbp = B.P.size();

        if (i == 0) {
            B.P.resize(nbp + np);

            for (j = 0; j < np; j++)
                B.P[j+nbp] = wall[i].P[j];
        }

        else {
            B.P.resize(nbp + np-1);

            for (j = 1; j < np; j++) 
                B.P[j-1+nbp] = wall[i].P[j];
        }  
    }
}

void createClosedBoundary(Wall& B, vector<Wall>& wall) {
    int i, j, np, nbp, nw = wall.size();

    for (i = 0; i < nw; i++) {
        np  = wall[i].P.size();
        nbp = B.P.size();

        if (i == 0) {
            B.P.resize(nbp + np);

            for (j = 0; j < np; j++)
                B.P[j+nbp] = wall[i].P[j];
        }

        else if (i == nw-1) {
            B.P.resize(nbp + np-2);
            for (j = 1; j < np-1; j++)
                B.P[j-1+nbp] = wall[i].P[j];
		}

        else {
            B.P.resize(nbp + np-1);

            for (j = 1; j < np; j++) 
                B.P[j-1+nbp] = wall[i].P[j];
        }  
    }
}

void connectAdjacentBoundaries(Wall& B1, Wall& B2) {
	int i, icase;

	int np1 = B1.P.size();
    int np2 = B2.P.size();

	findConnectParticle(B1, B2, icase);

	if (icase == 0) {
		B1.P.resize(np1+np2-1);

		sumUnitVector(B1.P[0], B2.P[0]);

		for (i = 0; i < np2-1; i++)
			B1.P[np1+i] = B2.P[(np2-1)-i];

	}

	else if (icase == 1) {
		
		//cout << "x = " << B1.P[np1-1].x << "\t y = " << B1.P[np1-1].y << endl;
		//cout << "x = " << B2.P[np2-1].x << "\t y = " << B2.P[np2-1].y << endl;
		//cout << "nx = " << B1.P[np1-1].nx << "\t ny = " << B1.P[np1-1].ny << endl;
		//cout << "nx = " << B2.P[np2-1].nx << "\t ny = " << B2.P[np2-1].ny << endl;
		
		B1.P.resize(np1+np2-1);
		

		//cout << "po souctu" << endl;
		//cout << "x = " << B1.P[np1-1].x << "\t y = " << B1.P[np1-1].y << endl;
		//cout << "nx = " << B1.P[np1-1].nx << "\t ny = " << B1.P[np1-1].ny << endl;
		//cout << "x = " << B2.P[np2-1].x << "\t y = " << B2.P[np2-1].y << endl;
		//cout << "nx = " << B2.P[np2-1].nx << "\t ny = " << B2.P[np2-1].ny << endl;
		//cout << "icase = 1" << endl;
		//cout << endl;
		

		for (i = 1; i < np2; i++)
			B1.P[np1+i-1] = B2.P[(np2-1)-i];
			//B1.P[np1+i-1].x = B1.P[np1+i-1].x + B2.P[(np2-1)-i].x;
			//B1.P[np1+i-1].y = B1.P[np1+i-1].y + B2.P[(np2-1)-i].y;
		sumUnitVector(B1.P[np1-1], B2.P[np2-1]);
		//sumUnitVector(B1.P[np1], B2.P[np2-1]);
	}

	else if (icase == 2) {
		B1.P.resize(np1+np2-1);

		sumUnitVector(B1.P[0], B2.P[np2-1]);

		for (i = 0; i < np2-1; i++)
			B1.P[np1+i] = B2.P[i]; 

		//sumUnitVector(B1.P[np1], B2.P[0]);
	}

	else if (icase == 3) {
		sumUnitVector(B1.P[np1-1], B2.P[0]);
		B1.P.resize(np1+np2-1);

		//cout << "icase = 3" << endl;

		for (i = 1; i < np2; i++)
			B1.P[np1+i-1] = B2.P[i];
		
		//sumUnitVector(B1.P[np1], B2.P[0]);
	}

	else if (icase == 4) {
		B1.P.resize(np1+np2-2);

		for (i = 1; i < np2-1; i++)
			B1.P[np1+i-1] = B2.P[(np2-1)-i]; 
		
		//sumUnitVector(B1.P[np1], B2.P[np2-1]);
		sumUnitVector(B1.P[0], B2.P[0]);
		sumUnitVector(B1.P[np1-1], B2.P[np2-1]);
	}

	else if (icase == 5) {
		B1.P.resize(np1+np2-2);

		for (i = 1; i < np2-1; i++)
			B1.P[np1+i-1] = B2.P[i]; 

		sumUnitVector(B1.P[np1-1], B2.P[0]);
		sumUnitVector(B1.P[0], B2.P[np2-1]);

	
		//sumUnitVector(B1.P[np1], B2.P[0]);
	}

}

void connectAllBoundaries(vector<Wall>& Bndr) {
	int n = Bndr.size();

	for (int i = 1; i < n; i++) 
		connectAdjacentBoundaries(Bndr[0], Bndr[i]);

}


void linspaceInterval(vector<double>& t, double t_0, double dt, double t_end, double& zb) {
	int n = (t_end - t_0)/dt + 1;
	zb = t_end - (t_0 + dt * (n-1));

	if (zb <= dt/2.) {
		n--;
		zb += dt;
	}

	t.resize(n);

	for (int i = 0; i < n; i++)
		t[i] = t_0 + i * dt; 

}

void linspaceInterval_fit(vector<double>& t, double t_0, double dt, double t_end, double r) {
	int n = (t_end - t_0)/dt + 1;
	double zb = t_end - (t_0 + dt * (n-1));

	double dtold = dt;
	double deltadt = zb/(n-1);
	dt += deltadt;

	cout << endl << "Note: wall nelze rozdelit rovnomerne s dp = " << dtold * r << endl;
	cout << "wall rozdelena rovnomerne s dp = " << dt * r << endl << endl;

	t.resize(n);

	for (int i = 0; i < n; i++)
		t[i] = t_0 + i * dt; 
}

bool compareFloatNumbers(double x, double y, double eps) {
	return fabs(x-y) <= eps;
}

void sumUnitVector(particle& P1, particle& P2) {
	double nx, ny, l;

	nx = P1.nx + P2.nx;
	ny = P1.ny + P2.ny;
	//cout << "nx po souctu = " << nx << endl;
	//cout << "ny po souctu = " << ny << endl;

	l  = sqrt(pow(nx, 2) + pow(ny, 2));
	if (l < 1e-15)
		l = 1;
	//cout << "l = " << l << endl;;

	P1.nx = P2.nx = nx / l;
	P1.ny = P2.ny = ny / l;
}

void findConnectParticle(Wall B1, Wall B2, int& icase) {
    double eps = 1e-14;

    int n1 = B1.P.size();
	int n2 = B2.P.size();

	if (compareFloatNumbers(B1.P[0].x, B2.P[0].x, eps) && compareFloatNumbers(B1.P[0].y, B2.P[0].y, eps) && 
			 compareFloatNumbers(B1.P[n1-1].x, B2.P[n2-1].x, eps) && compareFloatNumbers(B1.P[n1-1].y, B2.P[n2-1].y, eps)) {

		icase = 4;
    } 

	else if (compareFloatNumbers(B1.P[0].x, B2.P[n2-1].x, eps) && compareFloatNumbers(B1.P[0].y, B2.P[n2-1].y, eps) && 
			 compareFloatNumbers(B1.P[n1-1].x, B2.P[0].x, eps) && compareFloatNumbers(B1.P[n1-1].y, B2.P[0].y, eps)) {

		icase = 5;
    } 

    else if (compareFloatNumbers(B1.P[0].x, B2.P[0].x, eps) && 
		compareFloatNumbers(B1.P[0].y, B2.P[0].y, eps)) {
        
		icase = 0;
    }

    else if (compareFloatNumbers(B1.P[n1-1].x, B2.P[n2-1].x, eps) && 
			 compareFloatNumbers(B1.P[n1-1].y, B2.P[n2-1].y, eps)) {
        
		icase = 1;
    }

    else if (compareFloatNumbers(B1.P[0].x, B2.P[n2-1].x, eps) && 
			 compareFloatNumbers(B1.P[0].y, B2.P[n2-1].y, eps)) {

		icase = 2;
    }

    else if (compareFloatNumbers(B1.P[n1-1].x, B2.P[0].x, eps) && 
			 compareFloatNumbers(B1.P[n1-1].y, B2.P[0].y, eps)) {

		icase = 3;
    } 

	else
		cout << "Error: Bndr[" << B1.idx << "] a Bndr[" <<  B2.idx << "] nesdílí žádné společné částice" << endl;  
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

double integrateFunction(double f(double x), double a, double b) {
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