#include <vector>
#include "wall.h"


void define_normals(std::vector<Wall>& wall);

int main(){
	double dp = 1;

	std::vector<Wall> wall = {Wall(2, 2, 2, 8, dp), 
						 	  Wall(2, 8, 8, 8, dp),
						 	  Wall(8, 8, 8, 2, dp),
						 	  Wall(8, 2, 2, 2, dp)};

	define_normals(wall);
	
	return 0;
}

void define_normals(std::vector<Wall>& wall) {
	for (int i = 0; i < wall.size(); i++)
		wall[i].Save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");
}


