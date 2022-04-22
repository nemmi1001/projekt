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

	Wall(double x_0, double y_0, double x_end, double y_end, double dp);

    void Save(std::vector<particle>& P, const std::string& filename);
};

void Corners(std::vector<particle>& P1, std::vector<particle>& P2);
bool CompareFloatNumbers(double x, double y, double eps);

