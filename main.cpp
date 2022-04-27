#include <vector>
#include <iostream>
#include "wall.h"


double f(double x);


int main(){
	double dp = 0.5;

	std::vector<Wall> wall;
	//defineRectangle(wall, 2, 2, 6, 3, dp);
	//defineCircle(wall, 3, 5, 2, dp);
	defineCircleArc(wall, 6, 5, 1, 0, 6, 0, dp);
/*
	std::vector<Wall> wall = {Wall(2, 2, 2, 5, dp), 
						 	  Wall(2, 5, 6, 8, dp),
						 	  Wall(6, 8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};

	std::vector<Wall> wall = {Wall(2, 2, 2, 5.5, dp), 
						 	  Wall(2, 5.5, 6.3, 7.8, dp),
						 	  Wall(6.3, 7.8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};
*/
	
/*	
	int nsten = 6;
	wall.resize(nsten);
	
	wall[0].create(4, 2, 4, 3, dp);
	wall[1].create(4, 3, 2, 3, dp);
	wall[2].create(2, 3, 2, 6, dp);
	wall[3].create(2, 6, 6, 6, dp);
	wall[4].create(6, 6, 6, 2, dp);
	wall[5].create(6, 2, 4, 2, dp);
	*/

	//solveClosedWallCorners(wall);

	WallFinalize(wall);

	WallSave(wall);

	return 0;
}


double f(double x){
	return x*x;
}
