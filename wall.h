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
    int n;        // pocet castic steny
    double l;     // delka steny
    double s[2];  // smerovy vektor steny

    std::vector<particle> P;    // castice steny

	Wall(double x_0, double y_0, double x_end, double y_end, double dp);

    void Save(std::vector<particle>& P, const std::string& filename);
};

