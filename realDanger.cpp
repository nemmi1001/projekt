void prekvapko(std::vector<Wall>& Bndr, double dp) {
    Bndr.resize(7);

	std::vector<Wall> wall;
	std::vector<Wall> wall2;
	std::vector<Wall> wall3;
	std::vector<Wall> wall4(1);
	std::vector<Wall> wall5;
	std::vector<Wall> wall6(1);
	std::vector<Wall> wall7;

	defineCircleArcAlt(wall, 3, 3, 5, 3, 1, 1, 1, dp);    
    defineCircleArcAlt(wall2, 5, 3, 7, 3, 1, 1, 1, dp);	
    defineCircleArcAlt(wall3, 7, 3, 6, 4, 1, -1, 1, dp);

	wall4[0].create(6, 7, 6, 4, dp);

	defineCircleArcAlt(wall5, 6, 7, 4, 7, 1, 1, 1, dp);

	wall6[0].create(4, 4, 4, 7, dp);

	defineCircleArcAlt(wall7, 4, 4, 3, 3, 1, -1, 1, dp);
	
	createPartBoundary(Bndr[0], wall);
	createPartBoundary(Bndr[1], wall2);
	createPartBoundary(Bndr[2], wall3);
	createPartBoundary(Bndr[3], wall4);
	createPartBoundary(Bndr[4], wall5);
	createPartBoundary(Bndr[5], wall6);
	createPartBoundary(Bndr[6], wall7);

	//saveEachBoundary(Bndr);

	connectAllBoundaries(Bndr);
}