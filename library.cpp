#include "library.h"


void tvar1(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(1);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3(1);
	std::vector<Wall> wall4;

	wall[0].create(8, 2, 5, 2, dp);

	defineCircleArcAlt(wall2, 5, 2, 2, 5, 3, -1, -1, dp);  

    wall3[0].create(2, 5, 2, 8, dp);
    
	defineCircleArcAlt(wall4, 8, 2, 2, 8, 6, -1, 1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar2(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(3);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3(3);
	std::vector<Wall> wall4;

	wall[0].create(3, 6, 2, 6, dp);
	wall[1].create(2, 6, 2, 8, dp);
	wall[2].create(2, 8, 3, 8, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 8, 3, 3, 8, 5, -1, 1, dp);  

    wall3[0].create(8, 3, 8, 2, dp);
	wall3[1].create(8, 2, 6, 2, dp);
	wall3[2].create(6, 2, 6, 3, dp);
    solveNWallCorners(wall3);

	defineCircleArcAlt(wall4, 6, 3, 3, 6, 3, -1, -1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar3(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(2);

	std::vector<Wall> wall(4);
	std::vector<Wall> wall2;

	wall[0].create(7, 4, 7, 2, dp);
	wall[1].create(7, 2, 2, 2, dp);
	wall[2].create(2, 2, 2, 7, dp);
	wall[3].create(2, 7, 4, 7, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 4, 7, 7, 4, 3, 1, -1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar4(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(2);

	std::vector<Wall> wall(4);
	std::vector<Wall> wall2;

	wall[0].create(7, 4, 7, 2, dp);
	wall[1].create(7, 2, 2, 2, dp);
	wall[2].create(2, 2, 2, 7, dp);
	wall[3].create(2, 7, 4, 7, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 7, 4, 4, 7, 3, -1, 1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar5(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(1);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3(1);
	std::vector<Wall> wall4;

	wall[0].create(3, 2, 2, 2, dp);

	defineCircleArcAlt(wall2, 8, 2, 2, 2, 3, -1, 1, dp);

	wall3[0].create(8, 2, 7, 2, dp);

	defineCircleArcAlt(wall4, 7, 2, 3, 2, 2, -1, -1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar6(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(3);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3;
	std::vector<Wall> wall4;

	wall[0].create(8, 4, 8, 2, dp);
	wall[1].create(8, 2, 2, 2, dp);
	wall[2].create(2, 2, 2, 4, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 2, 4, 4, 6, 2, 1, -1, dp);
	defineCircleArcAlt(wall3, 6, 6, 4, 6, 1, -1, 1, dp);
	defineCircleArcAlt(wall4, 6, 6, 8, 4, 2, 1, -1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar7(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(3);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3(3);
	std::vector<Wall> wall4;

	wall[0].create(6, 7, 8, 7, dp);
	wall[1].create(8, 7, 8, 3, dp);
	wall[2].create(8, 3, 6, 3, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 6, 3, 4, 3, 2, -1, -1, dp);

	wall3[0].create(4, 3, 2, 3, dp);
	wall3[1].create(2, 3, 2, 7, dp);
	wall3[2].create(2, 7, 4, 7, dp);
	solveNWallCorners(wall3);

	defineCircleArcAlt(wall4, 4, 7, 6, 7, 2, 1, -1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar8(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(4);

	std::vector<Wall> wall(3);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3(3);
	std::vector<Wall> wall4;

	wall[0].create(6, 7, 8, 7, dp);
	wall[1].create(8, 7, 8, 3, dp);
	wall[2].create(8, 3, 6, 3, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 4, 3, 6, 3, 2, 1, 1, dp);

	wall3[0].create(4, 3, 2, 3, dp);
	wall3[1].create(2, 3, 2, 7, dp);
	wall3[2].create(2, 7, 4, 7, dp);
	solveNWallCorners(wall3);

	defineCircleArcAlt(wall4, 6, 7, 4, 7, 2, -1, 1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar9(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(1);

	std::vector<Wall> wall(4);

	wall[0].create(3, 2, 2, 7, dp);
	wall[1].create(2, 7, 8, 7, dp);
	wall[2].create(8, 7, 7, 2, dp);
	wall[3].create(7, 2, 3, 2, dp);
	solveClosedWallCorners(wall);

	createPartBoundary(Bndr[0], wall);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar10(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(1);

	std::vector<Wall> wall(4);

	wall[0].create(3, 2, 2, 4, dp);
	wall[1].create(2, 4, 8, 7, dp);
	wall[2].create(8, 7, 6, 2, dp);
	wall[3].create(6, 2, 3, 2, dp);
	solveClosedWallCorners(wall);

	createPartBoundary(Bndr[0], wall);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar11(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(1);

	std::vector<Wall> wall(6);

	wall[0].create(1, 2, 4, 6, dp);
	wall[1].create(4, 6, 4, 8, dp);
	wall[2].create(4, 8, 6, 8, dp);
	wall[3].create(6, 8, 6, 6, dp);
	wall[4].create(6, 6, 9, 2, dp);
	wall[5].create(9, 2, 1, 2, dp);
	solveClosedWallCorners(wall);

	createPartBoundary(Bndr[0], wall);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}

void tvar12(std::vector<Wall>& Bndr, double dp) {
	Bndr.resize(3);

	std::vector<Wall> wall(2);
	std::vector<Wall> wall2;
	std::vector<Wall> wall3;

	wall[0].create(8, 6, 5, 2, dp);
	wall[1].create(5, 2, 2, 6, dp);
	solveNWallCorners(wall);

	defineCircleArcAlt(wall2, 5, 6, 2, 6, 1.5, -1, 1, dp);
	defineCircleArcAlt(wall3, 8, 6, 5, 6, 1.5, -1, 1, dp);

	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);


	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}