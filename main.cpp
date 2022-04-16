#include "wall.h"

void define_normals(Wall wall[], const int size);

int main(){
	double dp = 1;

	const int n_wall = 4;
	Wall wall[n_wall] = {Wall(2, 2, 2, 8, dp), 
						 Wall(2, 8, 8, 8, dp),
						 Wall(8, 8, 8, 2, dp),
						 Wall(8, 2, 2, 2, dp)};

	for (int i = 0; i < n_wall; i++)
		wall[i].Save(wall[i].P, "wall" + std::to_string(i + 1) + ".txt");

	define_normals(wall, n_wall);
	
	return 0;
}

void define_normals(Wall wall[], const int size) {

}


