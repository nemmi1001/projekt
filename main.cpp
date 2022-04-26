#include <vector>
#include <iostream>
#include "wall.h"


double f(double x);


int main(){
	double dp = 0.1;

	std::vector<Wall> wall;
	//defineRectangle(wall, 2, 2, 6, 3, dp);
       //	defineCircle(wall, 3, 5, 2, dp);
       defineCircleArc(wall,5,5,-5,5,0,0,dp);

/*
	std::vector<Wall> wall = {Wall(2, 2, 2, 5, dp), 
						 	  Wall(2, 5, 6, 8, dp),
						 	  Wall(6, 8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};
	solveCorners(wall);
*/
	WallFinalize(wall);

	WallSave(wall);

	return 0;
}


double f(double x){
	return x*x;
}
