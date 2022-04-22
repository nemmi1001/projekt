#include <vector>
#include "wall.h"


void define_normals(std::vector<Wall>& wall);

int main(){
	double dp = 1;

	std::vector<Wall> wall;
	defineRectangle(wall, 2, 2, 6, 6, dp);
	//defineCircle(wall, 5, 5, 3, dp);
/*
	std::vector<Wall> wall = {Wall(2, 2, 2, 5, dp), 
						 	  Wall(2, 5, 6, 8, dp),
						 	  Wall(6, 8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};
	solveCorners(wall);
*/
	saveMesh(wall);

	return 0;
}




