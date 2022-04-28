#include <vector>
#include <iostream>
#include "wall.h"

double f(double x);


int main(){
	double dp = 0.1;
	//std::vector<Wall> wall;
	//defineRectangle(wall, 2, 2, 6, 3, dp);
       //	defineCircle(wall, 3, 5, 2, dp);
       //defineCircleArc(wall,5,5,-5,5,0,0,dp);
	//defineCADCircleArc(wall,1,0,3,0,1,-1,dp);

	std::vector<Wall> wall = {Wall(2, 2, 2, 5, dp), 
						 	  Wall(2, 5, 6, 8, dp),
						 	  Wall(6, 8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};
	solveWallCorners(wall);

	WallFinalize(wall);

	WallSave(wall);
   	write_to_ASCII_VTK(wall[0],"gg.vtk");
	return 0;
}


double f(double x){
	return x*x;
}
