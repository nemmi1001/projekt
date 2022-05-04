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

    void create(double x_0, double y_0, double x_end, double y_end, double dp);         // doplni castici, ktera nebude vzdalena o dp
    void create_fit(double x_0, double y_0, double x_end, double y_end, double dp);     // upravi dp, aby byly castice vzdaleny rovnomerne
    
    void save(std::vector<particle>& P, const std::string filename);
    void save2VTK(std::vector<particle>& P, const std::string filename);
};

void solveWallCorner(Wall& w1, Wall& w2);
void solveNWallCorners(std::vector<Wall>& wall);
void solveClosedWallCorners(std::vector<Wall>& wall);

void defineRectangle(std::vector<Wall>& wall, double x_0, double y_0, double a, double b, double dp);
void defineCircle(std::vector<Wall>& wall, double x_0, double y_0, double r, double dp);
void defineCircle_fit(std::vector<Wall>& wall, double x_0, double y_0, double r, double dp);
void defineCircleArc(std::vector<Wall>& wall, double x_0, double y_0, double x_end, double y_end, double x_s, double y_s, int ipar, double dp);
void defineCircleArcAlt(std::vector<Wall>& wall, double x_0, double y_0, double x_b, double y_b, double r, int ori, int ipar, double dp);

void createPartBoundary(Wall& B, std::vector<Wall>& wall);
void createClosedBoundary(Wall& B, std::vector<Wall>& wall);
void connectAdjacentBoundaries(Wall& B1, Wall& B2);
void connectAllBoundaries(std::vector<Wall>& Bndr);

void saveEachWall(std::vector<Wall>& wall);
void saveEachBoundary(std::vector<Wall>& Bndr);

bool compareFloatNumbers(double x, double y, double eps);
void linspaceInterval(std::vector<double>& t, double t_0, double dt, double t_end, double& zb);
void linspaceInterval_fit(std::vector<double>& t, double t_0, double dt, double t_end, double r);
void sumUnitVector(particle& P1, particle& P2);
void findConnectParticle(Wall B1, Wall B2, int& icase);

double integrateRectMethod(double x, double dx, double f(double x));
double integrateTrapMethod(double x, double dx, double f(double x), double a, double b);
double integrateFunction(double f(double x), double a, double b);
