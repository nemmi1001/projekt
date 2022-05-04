#include <vector>
#include <iostream>
#include "wall.h"
#include "library.h"


int main(){
	double dp = 0.25;

/* ---------- DEFINOVANÉ UZAVŘENÉ TVARY --------- *//*
	Wall Bndr;
	std::vector<Wall> wall;

	//defineRectangle(wall, 2, 2, 6, 3, dp); 		// ([x0, y0], a, b)
	defineCircle(wall, 5, 5, 2, dp);			// ([xS, yS], r)

	createClosedBoundary(Bndr, wall);

	Bndr.save(Bndr.P, "closed_boundary.txt");
	Bndr.save2VTK(Bndr.P, "closed_boundary.vtk");


*//* ---------------- VLASTNÍ TVARY --------------- */
	/* ----------------- ROVINNÉ STĚNY ------------------- */
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


	/* -------- KOMBINACE ROVINNÉ STĚNY A KRUHU -------- *//*
	std::vector<Wall> Bndr(4);

	std::vector<Wall> wall;
	std::vector<Wall> wall2(2);
	std::vector<Wall> wall3;
	std::vector<Wall> wall4(2);

	defineCircleArc(wall, 5, 3, 7, 3, 6, 3, 1, dp); 		// ([x0, y0], [xend, yend], [xS, yS], normála = 1 / -1)
	//defineCircleArcAlt(wall, 5, 3, 7, 3, 1, -1, 1, dp);	// ([x0, y0], [xend, yend], r, orientace = 1 / -1, normála = 1 / -1)
	
	wall2[0].create(3, 7, 7, 7, dp);
	wall2[1].create(7, 7, 7, 3, dp);
	solveNWallCorners(wall2);
	
	defineCircleArc(wall3, 3, 5, 3, 7, 3, 6, -1, dp);
	//defineCircleArcAlt(wall3, 3, 5, 3, 7, 1, 1, -1, dp);
	
	wall4[0].create(5, 3, 5, 5, dp);
	wall4[1].create(5, 5, 3, 5, dp);
	solveNWallCorners(wall4);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);

	Bndr[0].save(Bndr[0].P, "closed_boundary.txt");
	Bndr[0].save2VTK(Bndr[0].P, "closed_boundary.vtk");


	*//* ---------- KONKRÉTNÍ DEFINOVANÉ TVARY --------------- *//*
	std::vector<Wall> Bndr;

	tvar1(Bndr, dp);
	//tvar2(Bndr, dp);
	//tvar3(Bndr, dp);
	//tvar4(Bndr, dp);
	//tvar5(Bndr, dp);
	//tvar6(Bndr, dp);
	//tvar7(Bndr, dp);
	//tvar8(Bndr, dp);
	//tvar9(Bndr, dp);
	//tvar10(Bndr, dp);
	//tvar11(Bndr, dp);
	//tvar12(Bndr, dp);

	Bndr[0].save(Bndr[0].P, "closed_boundary.txt");
	Bndr[0].save2VTK(Bndr[0].P, "closed_boundary.vtk");
	*//* ---------------------------------------------------- */

	return 0;
}