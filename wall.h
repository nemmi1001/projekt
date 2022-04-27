#pragma once

#include <vector>
#include <string>

typedef struct{
	double x, y, nx, ny;
} particle;


class Wall {
private:
	/* data */
public:
    int idx;      // index steny
    int n;        // pocet castic steny
    double l;     // delka steny
    double s[2];  // smerovy vektor steny

    std::vector<particle> P;    // castice steny

    Wall(){};
	Wall(double x_0, double y_0, double x_end, double y_end, double dp);

    void create(double x_0, double y_0, double x_end, double y_end, double dp);
    void create_fit(double x_0, double y_0, double x_end, double y_end, double dp);     // upravi dp, aby byly castice vzdaleny rovnomerne
    void create_insert(double x_0, double y_0, double x_end, double y_end, double dp);  // doplni castici, ktera nebude vzdalena o dp
    void save(std::vector<particle>& P, const std::string& filename);
};

void solveClosedWallCorners(std::vector<Wall>& wall);
void solveWallCorner(Wall& w1, Wall& w2);

void defineRectangle(std::vector<Wall>& wall, double x_0, double y_0, double a, double b, double dp);
void defineCircle(std::vector<Wall>& wall, double x_0, double y_0, double r, double dp);
void defineCircleArc(std::vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double x_s, double y_s, double dp);

void WallSave(std::vector<Wall>& wall);
void WallFinalize(std::vector<Wall>& wall);


bool compareFloatNumbers(double x, double y, double eps);
void linspaceInterval(std::vector<double>& t, double t_0, double dt, double t_end, double* zb);
void sumUnitVector(particle& P1, particle& P2);
double integrateRectMethod(double x, double dx, double f(double x));
double integrateTrapMethod(double x, double dx, double f(double x), double a, double b);
double integrateFunction(double f(double x), double a, double b);
