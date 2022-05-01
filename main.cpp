#include <vector>
#include <iostream>
#include "wall.h"


double f(double x);


int main(){
	double dp = 0.5;

/* ---------- DEFINOVANÉ UZAVŘENÉ TVARY --------- */
	Wall Bndr;
	std::vector<Wall> wall;

	defineRectangle(wall, 2, 2, 6, 3, dp); 		// ([x0, y0], a, b)
	//defineCircle(wall, 3, 5, 2, dp);			// ([xS, yS], r)
	
	createClosedBoundary(Bndr, wall);

	Bndr.save(Bndr.P, "closed_boundary.txt");
	Bndr.save2VTK(Bndr.P, "closed_boundary.vtk");


/* ---------------- VLASTNÍ TVARY --------------- */
	/* ----------------- ROVINNÉ STĚNY ------------------- *//*
	Wall Bndr;
	std::vector<Wall> wall(6);

	wall[0].create(4, 2, 4, 3, dp); 	// ([x0, y0], [xend, yend])
	wall[1].create(4, 3, 2, 3, dp);
	wall[2].create(2, 3, 2, 6, dp);
	wall[3].create(2, 6, 6, 6, dp);
	wall[4].create(6, 6, 6, 2, dp);
	wall[5].create(6, 2, 4, 2, dp);	

	solveClosedWallCorners(wall);

	createClosedBoundary(Bndr, wall);

	Bndr.save(Bndr.P, "closed_boundary.txt");
	Bndr.save2VTK(Bndr.P, "closed_boundary.vtk");


	*//* ----- KOMBINACE ROVINNÉ STĚNY A KRUHOVÉ ČÁSTI ----- *//*
	std::vector<Wall> Bndr(4);

	std::vector<Wall> wall(6);
	//defineCircleArc(wall, 7, 6, 2, 1, 7, 1, dp); 		// ([x0, y0], [xend, yend], [xS, yS])
	defineCircleArcAlt(wall, 5, 3, 7, 3, 1, -1, dp);	// ([x0, y0], [xend, yend], r, orientace = 1, -1)	
	createPartBoundary(Bndr[0], wall);

	std::vector<Wall> wall2(2);
	wall2[0].create(3, 7, 7, 7, dp);
	wall2[1].create(7, 7, 7, 3, dp);
	solveNWallCorners(wall2);
	createPartBoundary(Bndr[1], wall2);

	std::vector<Wall> wall3;
	defineCircleArcAlt(wall3, 3, 7, 3, 5, 1, -1, dp);
	createPartBoundary(Bndr[2], wall3);

	std::vector<Wall> wall4(2);
	wall4[0].create(5, 3, 5, 5, dp);
	wall4[1].create(5, 5, 3, 5, dp);
	solveNWallCorners(wall4);
	createPartBoundary(Bndr[3], wall4);

	saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);

	Bndr[0].save(Bndr[0].P, "closed_boundary.txt");
	Bndr[0].save2VTK(Bndr[0].P, "closed_boundary.vtk");
	*//* --------------------------------------------------- */

	return 0;
}


double f(double x){
	return x*x;
}
