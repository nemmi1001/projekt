#include <vector>
#include "wall.h"


void define_normals(std::vector<Wall>& wall);

int main(){
	double dp = 1;

	std::vector<Wall> wall = {Wall(2, 2, 2, 5, dp), 
						 	  Wall(2, 5, 6, 8, dp),
						 	  Wall(6, 8, 6, 2, dp),
						 	  Wall(6, 2, 2, 2, dp)};
	Corners(wall[0].P, wall[1].P);
	Corners(wall[1].P, wall[2].P);
	Corners(wall[2].P, wall[3].P);
	Corners(wall[3].P, wall[0].P);

	define_normals(wall);
	
	return 0;
}

void define_normals(std::vector<Wall>& wall) {
	for (int i = 0; i < wall.size(); i++)
		wall[i].Save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");
}


