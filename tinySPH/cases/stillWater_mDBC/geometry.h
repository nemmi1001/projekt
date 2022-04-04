#include "SPH_defs.h"
#include "SPH_particle_system.h"

		// WIDE BOX
		Create_Wall_x_mdbc(particles, 0.+dp , 0.8- dp, 0. , -1, dp/2);
		Create_Wall_y_mdbc(particles, 0. , 0 + dp, .5 , -1, dp/2,1);
		Create_Wall_y_mdbc(particles, 0.8 , 0 + dp, .5 , +1, 0.8-dp/2,10000);
		Create_Wall_x_mdbc(particles, 0.+dp , 0.8- dp, .5+dp , 1, dp/2);
		
		Create_corner_mdbc(particles, -1, -1, {0. + dp/2, 0. + dp/2},1);
		Create_corner_mdbc(particles, +1, -1, {0.8 - dp/2, 0. + dp/2},10000);
		//moje
		Create_corner_mdbc(particles, -1, 1, {0. + dp/2, 0.5 + dp/2},1);
		Create_corner_mdbc(particles, 1, 1, {.8 - dp/2, 0.5 + dp/2},10000);
	        //Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.8, 0.5);
		//Create_solid(particles, 0 + dp, 0 + dp, 0.8-dp, 0.5);

		// NARROW BOX
		// Create_Wall_x_mdbc(particles, 0.+dp , 0.3, 0. , -1, dp/2);
		// Create_Wall_y_mdbc(particles, 0. , 0 + dp, 1. , -1, dp/2);
		// Create_Wall_y_mdbc(particles, 0.3 , 0 + dp, 1. , +1, 0.3-dp/2);

		// Create_corner_mdbc(particles, -1, -1, {0. + dp/2, 0. + dp/2});
		// Create_corner_mdbc(particles, +1, -1, {0.3 - dp/2, 0. + dp/2});

	        //Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		 Create_Block_of_FluidWP(particles, 0 + dp, 0 + dp, 0.8-dp, 0.5);
