#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
Ex 1
*/
#define N 100 // number of fluid zones to simulate
#define gamma 1.4 // adiabatic constant
void pdeint(); // main driver to integrate eulers eqs for a shock tube 
void boundary(); // handles updating of ghost cells 
void lax_fried_stepper(); // makes a step in time based on lax friedrichs finite volume scheme
void hhlc_rk_stepper(); // makes a step in time based on 3rd order rk and hhlc finite volume scheme
void main()
{	
	
	
	// allocate memory for state and flux arrays of shape (N,3), where N is the number of fluid zones.
	// each fluid zone has 3 elements refering to the density / density flux, momentum / momentum flux, energy / energy flux, respectively.
	double (*Q)[3] = malloc(sizeof(double[N][3])); 


	
// wrap this into func : args = [starting q and f, start and end times, saveflag, q and fvecs,integration method, derivs]
	//run sim
	t = tstart;

	if (saveflag)
	{
		char outfile[] = "log_1.csv"; // log file name
		FILE *pf; 
		pf = fopen(outfile,"w+");
	}
 
	while (t<tend)
	{
		//[STEP 1] compute fluxes within cells
		
		//[STEP 2] determine fluxes between cells

		//[STEP 3] update state vectors

		for (int i=1; i<N-1; i++) qnew[i] = q[i];
		{
			// update boundary vals
			qnew[0] = qnew[1];
			qnew[-1] = qnew[-2];
		}

		memcpy(q,qnew,N*sizeof(double));

		//[STEP 4] update time

		//[STEP 5] save data?

	

	

		}
	}
}
	

 
void shock_pde_int(state_array, t_start, t_end, x_start, x_end, N, boundary_pos, boundary_left, boundary_right, boundary_condition, gamma, saveflag, stepper)
{
	/* Integrates eulers equations to solve 1d shock tube problem given initial conditions
	and choice of stepper routine. state_array will store values of fluid cells at end time.
	

	Parameters:
	-----------
		state_array: 
			2d fluid cell state array of doubles with shape (n_fluid_cells,3), col headers = [rho, rho*v, e], where rho is fluid density, rho*v is momentum density, e is energy density. 
		
		t_start: double
			starting time of simulation
		
		t_end: double
			ending time of simulation

		x_start : double
			x position of first fluid cell

		x_end : double
			x position of last fluid cell
		
		N: int 
			number of fluid cells excluding ghost cells
		
		boundary_pos: double
			x position of shock tube boundary at start_time
		
		boundary_left: array
			array of boundary conditions for left side of shocktube. Array expected to contain entries [rho, p, v], where rho is fluid density, p is pressure, v is velocity.
		
		boundary_right: array
			array of boundary conditions for right side of shocktube. Array expected to contain entries [rho, p, v], where rho is fluid density, p is pressure, v is velocity.
		
		boundary_condition: str
			str specifying how boundaries / ghost cells are updated. Choices = ['outflow','periodic', 'reflecting']
		
		gamma: double
			value of fluids polytropic index 
		
		saveflag: int
			set to a non zero value to save integrated fluid cell states at each time step to log.txt file
		
		stepper: function
			external function which is called at every time step to update cell state and cell flux arrays.
	
	Returns:
	--------
		None
	*/
	/*TODO:
	understand differences betweeen lax friend, hll, hllc etc to figure out how to handle updates

	decide if data saved at current time or current time + dt per iteration
	*/

	// define independant and dependant variables along with step size vars
	double dx,dt,t,x; // grid_cell_size, time step, time and position vars;
	double v,rho,p,e; // velocity, density, pressure, energy density 
	dx = (x_end-x_start)/N; 
	t = tstart;
	
	// define (N,3) array to keep track of cell fluxes at each time
	double (*F)[3] = malloc(sizeof(double[N][3]));

	// unpack boundary arrays
	double rho_left = boundary_left[0];
	double p_left = boundary_left[1];
	double v_left = boundary_left[2];
	double e_left = p_left / (gamma - 1);

	double rho_right = boundary_right[0];
	double p_right = boundary_right[1];
	double v_right = boundary_right[2];
	double e_right = p_right / (gamma - 1);

 	// [STEP 0] set initial values for state and flux arrays: q and f
	for (int i=0; i<N; i++)
	{
		x = i*dx; 
		if (x<x_boundary)
		{
			// set initial left state
			Q[i][0] = rho_left;
			Q[i][1] = rho_left * v_left;
			Q[i][2] = e_left;
		} 
		
		else 
		{
			// initial right state
			Q[i][0] = rho_right;
			Q[i][1] = rho_right * v_right;
			Q[i][2] = e_right;
		}
		
 	}

	// run sim
	while (t<t_end)
	{
		// [STEP 1] determine time step dt
		dt = ;
		
		// [STEP 2] set fluxes within cells
		for (int i=0; i<N; i++)
		{	
			p  = (gamma - 1) * Q[i][2]
			F[i][0] = Q[i][1];
			F[i][1] = pow(Q[i][1],2)/Q[i][0] + p;
			F[i][2] = Q[i][1] * (Q[i][2] + p_left) / Q[i][0];

			// printf("cell %i: \ninitial state:(%lf,%lf,%lf)\ninitial fluxes:(%lf,%lf,%lf)\n\n",i,q[i][0],q[i][1],q[i][2],f[i][0],f[i][1],f[i][2]);
		}

		// [STEP 3] update cell state array Q using stepper
		stepper(Q,F);

		// [STEP 4] update time
		dt = calc_cdf_step(Q);
		t += dt;
	}	
	free(F);
	// Save data?
	if (saveflag)
	{
		// save x,rho,p,v,e values to csv file
		for (int i=0; i<N; i++)
		{
			x = i*dx;
			fprintf(pf,"%lf,%lf\n", x, Q[i][0], (gamma - 1)*Q[i][2], Q[i][1]/Q[i][0], Q[i][2]);
		}
	}
}
