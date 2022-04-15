#include "wall.h"

void define_normals(Wall w1, Wall w2, Wall w3, Wall w4);

int main(){
	double dp = 1;

	Wall wall_1(2					  , 2					  , 2			 , 8			, dp), 
		 wall_2(wall_1.P[wall_1.n-1].x, wall_1.P[wall_1.n-1].y, 8			 , 8			, dp),
		 wall_3(wall_2.P[wall_2.n-1].x, wall_2.P[wall_2.n-1].y, 8			 , 2			, dp),
		 wall_4(wall_3.P[wall_3.n-1].x, wall_3.P[wall_3.n-1].y, wall_1.P[0].x, wall_1.P[0].y, dp);

	wall_1.Save(wall_1.P, "wall1.txt");
	wall_2.Save(wall_2.P, "wall2.txt");
	wall_3.Save(wall_3.P, "wall3.txt");
	wall_4.Save(wall_4.P, "wall4.txt");

	define_normals(wall_1, wall_2, wall_3, wall_4);
	return 0;
}

void define_normals(Wall w1, Wall w2, Wall w3, Wall w4) {

}


