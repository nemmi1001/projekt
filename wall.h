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
    void save(std::vector<particle>& P, const std::string& filename);
};

void solveCorners(std::vector<Wall>& wall);
void Corner(std::vector<particle>& P1, std::vector<particle>& P2);
bool compareFloatNumbers(double x, double y, double eps);
void defineNormals(std::vector<Wall>& wall);
void defineRectangle(std::vector<Wall>& wall, double x_0, double y_0, double a, double b, double dp);
void defineCircle(std::vector<Wall>& wall, double x_0, double y_0, double r, double dp);
void saveMesh(std::vector<Wall>& wall);

