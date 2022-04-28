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

void defineCircleArc(vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double x_s, double y_s, double dp) {
	vector <double> t;	
	double r=sqrt((x_0-x_s)*(x_0-x_s)+(y_0-y_s)*(y_0-y_s));
	double l_half=0.5*sqrt((-x_0+x_end)*(-x_0+x_end)+(-y_0+y_end)*(-y_0+y_end));
	double tau=2*asin(l_half/r);
	double psi;
	if(x_s==x_0 && y_0>y_s){psi=0.5*M_PI;} //horní pi/2
	else if (x_s==x_0 && y_0<y_s){psi=1.5*M_PI;}//dolní pi/2
	else {
	double psi_p=atan(fabs(y_0-y_s)/fabs(x_0-x_s));
	
	if(x_s<x_0 && y_s<=y_0){psi=psi_p;}//1. kvadrant
	else if(x_s>=x_0 && y_s<y_0){psi=M_PI-psi_p;}//2.kvadrant
	else if(x_s>x_0 && y_s>=y_0){psi=M_PI+psi_p;}//3.kvadrant
	else if(x_s<=x_0 && y_s>y_0){psi=2*M_PI-psi_p;}//4.kvadrant
	else {cout<<"nefunguje"<<endl;}
	
		}
	
	double phi =psi+tau;	
	//cout<<psi<< "   "<< phi <<endl;
	linspace(t,psi,dp/r,phi);
	int n= t.size();
	
	
	wall.resize(n);
		for(int i = 0; i < n-1; i++){
		//cout<<t[i]<<endl;
		wall[i].P.resize(2);
		wall[i].P[0].x  = x_s +r*cos(t[i]); 
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
	wall[n-1].P[1].nx =cos(phi);
	wall[n-1].P[1].ny =sin(phi);	
		
 


}


void defineCADCircleArc(vector<Wall>& wall, double x_0, double y_0, double x_b, double y_b, double r, int ori,double dp){
			double x_a = 0.5*(x_b+x_0);
			double y_a = 0.5*(y_0+y_b);  // průsečík osy úsečky s úsečkou
			double l_half = sqrt((x_a-x_0)*(x_a-x_0)+(y_a-y_0)*(y_a-y_0));
			double h = sqrt(r*r-l_half * l_half);
			
			//
			double x_1 = (y_a*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) - y_0*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) + h*h*x_0 - h*h*x_a - r*r*x_0 + r*r*x_a - x_0*x_a*x_a - x_0*x_0*x_a + x_0*y_0*y_0 + x_0*y_a*y_a + x_a*y_0*y_0 + x_a*y_a*y_a + x_0*x_0*x_0 + x_a*x_a*x_a - 2*x_0*y_0*y_a - 2*x_a*y_0*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double x_2 = (y_0*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) - y_a*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) + h*h*x_0 - h*h*x_a - r*r*x_0 + r*r*x_a - x_0*x_a*x_a - x_0*x_0*x_a + x_0*y_0*y_0 + x_0*y_a*y_a + x_a*y_0*y_0 + x_a*y_a*y_a + x_0*x_0*x_0 + x_a*x_a*x_a - 2*x_0*y_0*y_a - 2*x_a*y_0*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double y_1 = (x_0*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) - x_a*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) + h*h*y_0 - h*h*y_a - r*r*y_0 + r*r*y_a + x_0*x_0*y_0 + x_0*x_0*y_a + x_a*x_a*y_0 + x_a*x_a*y_a - y_0*y_a*y_a - y_0*y_0*y_a + y_0*y_0*y_0 + y_a*y_a*y_a - 2*x_0*x_a*y_0 - 2*x_0*x_a*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));

double y_2 = (x_a*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) - x_0*sqrt((h*h + 2*h*r + r*r - x_0*x_0 + 2*x_0*x_a - x_a*x_a - y_0*y_0 + 2*y_0*y_a - y_a*y_a)*(- h*h + 2*h*r - r*r + x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a)) + h*h*y_0 - h*h*y_a - r*r*y_0 + r*r*y_a + x_0*x_0*y_0 + x_0*x_0*y_a + x_a*x_a*y_0 + x_a*x_a*y_a - y_0*y_a*y_a - y_0*y_0*y_a + y_0*y_0*y_0 + y_a*y_a*y_a - 2*x_0*x_a*y_0 - 2*x_0*x_a*y_a)/(2*(x_0*x_0 - 2*x_0*x_a + x_a*x_a + y_0*y_0 - 2*y_0*y_a + y_a*y_a));
			//
			cout<<x_1<< "   "<<y_1<<"\n";
			cout<<x_2<< "   "<<y_2<<"\n";
			// ori=1  konvex, ori=-1 konkáv
			double px_konvex=0;
			double py_konvex=0;
			double px_konkav=0;
			double py_konkav=0;
			double p_x=0;
			double p_y=0;
			if (y_1>y_2){
				 py_konvex=y_1;
				 px_konvex=x_1;
				
				 py_konkav=y_2;
				 px_konkav=x_2;			}
			else if (y_1<=y_2){
				 py_konvex=y_2;
				 px_konvex=x_2;
				
				 py_konkav=y_1;
				 px_konkav=x_1;	}
			
			if(ori==1){
			 	 p_x=px_konvex;
				 p_y=py_konvex;
				
					;}
			else if(ori==-1){
				 p_x=px_konkav;
				 p_y=py_konkav;
					;}
			
			else {cout<<"Špatná orientace zadej 1 nebo -1"<<"\n";}
		
		
			
		
			defineCircleArc(wall, x_0, y_0, x_b , y_b,  p_x,  p_y, dp);
			





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
	//t[0] = t_0;

	for (int i = 0; i < n; i++)
		t[i] = t_0 + i * dt; 

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


//hokus pokus
void write_to_ASCII_VTK (Wall wall, string filename)
{

        int np= wall.P.size();
        

        ofstream vystup(filename);
        

        // file header
        vystup << "# vtk DataFile Version 3.0" << endl << "vtk output" << endl << "ASCII" << endl << "DATASET POLYDATA" << endl;

        // point coordinates list
        vystup << "POINTS " << np << " float" << std::endl;
        for(int i=0;i<np;i++)
        {
                vystup << wall.P[i].x << " " << wall.P[i].y << " 0" << endl;
        }
        
        // nefungují normály v paraview
      	vystup << "POLYDATA " << std::endl;

	/*vystup << "type 1 " << np << " int" << std::endl;
	for(int i=0;i<np;i++)
	{
		vystup << "f" << std::endl;
	}
	*/
	vystup << "norma " << np << " float" << std::endl;
	for(int i=0;i<np;i++)
	{
		vystup << wall.P[i].nx << wall.P[i].nx<< "0"<< std::endl;
	}
	
	

        //      data fields
        /*file << "POINT_DATA " << np << std::endl << "FIELDS FieldData 4" << std::endl;

        file << "type 1 " << np << " int" << std::endl;
        for(const idx &type: particles.data.part_type)
        {
                file << type << std::endl;
        }

        file << "pressure 1 " << np << " float" << std::endl;
        for(const real &p: particles.data.p)
        {
                file << p << std::endl;
        }

        file << "density 1 " << np << " float" << std::endl;
        for(const real &rho: particles.data.rho)
        {
                file << rho << std::endl;
        }

        file << "velocity 3 " << np << " float" << std::endl;
        for(const realvec &v : particles.data.v){
                file << v.x << " " << v.y << " 0" << std::endl;
        }
*/
        vystup.close();

}







