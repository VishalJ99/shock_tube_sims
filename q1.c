#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
Ex 1
*/

/*TODO:
write bc handler
Write hll
write hllc

Do flux arrays have ghost cells?
*/

#define N 100 // number of fluid zones to simulate
void shock_pde_int(double [N+2][3], double , double , double , double , int , double , double [3], double [3], int , double , int , void (*)(double [N+2][3], double [N+2][3], int, double, double)); // main driver to integrate eulers eqs for a shock tube 
void set_ghosts(double [N+2][3], int); // handles updating of ghost cells 
void lax_fried_stepper(double [N+2][3], double [N+2][3], int , double , double);
; // makes a step in time based on lax friedrichs finite volume scheme
void hhlc_rk_stepper(); // makes a step in time based on 3rd order rk and hhlc finite volume scheme
double get_cfl_dt(double [N+2][3], double , int , double , double );

// outfile
FILE *pf;

void main()
{	
	// memcpy(q,qnew,N*sizeof(double));
	double (*Q)[3] = malloc(sizeof(double[N+2][3]));
	double boundary_left[3] = {1, 1, 0.75};
	double boundary_right[3] = {0.125, 0.1, 0};
	double boundary_pos = 0.3;
	double t_start = 0;
	double t_end = 0.2;
	double x_start = 0;
	double x_end = 1;
	int saveflag = 1;
	double gamma = 1.4;
	int boundary_condition = 0;

	char outfile[] = "log.csv";
	pf = fopen(outfile,"w+");

	shock_pde_int(Q, t_start,t_end, x_start, x_end, N, boundary_pos, boundary_left, boundary_right, boundary_condition, gamma, saveflag, lax_fried_stepper);
}
void set_ghosts(double Q[N+2][3], int boundary_condition)
{
	if (boundary_condition == 0)
	{
		// outflow
		for (int j = 0; j<3; j++)
		{
			Q[0][j] = Q[1][j];
			Q[N+1][j] = Q[N][j];
		}
	}
	if (boundary_condition == 1)
	{
		// periodic 
		for (int j = 0; j<3; j++)
		{
			Q[0][j] = Q[N][j];
			Q[N+1][j] = Q[1][j];
		}
	}

	if (boundary_condition == 2)
	{
		// reflective 
		for (int j = 0; j<3; j++)
		{
			int sgn;
			sgn = (j==1) ? -1 : 1;
			
				Q[0][j] = sgn*Q[1][j];
				Q[N+1][j] = sgn*Q[N][j];
		}
	}
}	

double get_cfl_dt(double Q[N+2][3], double dx, int n_cells, double gamma, double margin)
{
	double dt,a,v,p,c,vmax,cmax;
	vmax = 0;
	cmax = 0;
	
	for (int i=0;i<n_cells+2;i++)
	{
		v = Q[i][1] / Q[i][0];
		p = (gamma - 1)*Q[i][2];
		c = pow(gamma*p/Q[i][0],0.5);
		if (abs(v) > vmax) vmax = v;
		if (c > cmax) cmax = c;
	}

	a = vmax + cmax; 
	dt = dx/a * margin;
	return dt;

}
void lax_fried_stepper(double Q[N+2][3], double F[N+2][3], int n_cells, double dx, double dt)
{
	/*
	Updates fluid cell state array Q using the Lax-Friedrichs finite volume scheme

	Parameters:
	-----------
		Q : array
			Fluid cell state array of shape (n_fluid_cells, 3). Each row contains entries [rho, rho*v, e] where rho is fluid density, v is velocity, e is internal energy.
		
		F : array
			Flux cell state array of shape (n_fluid_cells, 3). Each row contains entries
			of density flux, momentum_density flux, energy density flux within each cell i.
		
		n_cells: int
			Number of fluid cells

		dx : double
			size of each fluid cell

		dt : double
			step in time

	Returns:
	--------
		None
	*/

	// CHECK GHOST CELLS HANDLED CORRECTLY
	// CHECK DX DT VALUES CORRECT
	// CHECK F LEFT AND F RIGHT CORRECTY INDEXED AND SET

	// initialise left and right flux arrays
	double (*F_l)[3] = malloc(sizeof(double[N][3]));
	double (*F_r)[3] = malloc(sizeof(double[N][3]));
	// compute fluxes between cells (how to handle boundaries?) 
	for (int i = 0; i<n_cells; i++)
	{	
		for (int j = 0;  j<3; j++)
		{
			F_r[i][j] = 0.5 * ( (F[i][j] + F[i+1][j]) + (dx/dt) * (Q[i][j] - Q[i+1][j]) );  
			F_l[i][j] = 0.5 * ( (F[i-1][j] + F[i][j]) + (dx/dt) * (Q[i-1][j] - Q[i][j]) );  
		}
	}

	// update cell states
	for (int i = 1; i<=n_cells; i++)
	{
		for (int j = 0;  j<3; j++)
		{	
			// printf("before: Q[i][j] %lf, F_r[i][j] %lf, F_l[i][j] %lf\n",Q[i][j],F_r[i][j],F_l[i][j]);
			Q[i][j] = Q[i][j] - dt/dx * (F_r[i][j] - F_l[i][j]); 
			// printf("Q[I][J] after = %lf\n",Q[i][j]);
		}
	}
}
 
void shock_pde_int(double Q[N+2][3], double t_start, double t_end, double x_start, double x_end, int n_cells, double boundary_pos, double boundary_left[3], double boundary_right[3], int boundary_condition, double gamma, int saveflag, void (*stepper)(double [N+2][3], double [N+2][3], int, double, double))
{
	/* Integrates eulers equations to solve 1d shock tube problem given initial conditions
	and choice of stepper routine. state_array will store values of fluid cells at end time.
	

	Parameters:
	-----------
		Q: double array
			2d fluid cell state array of doubles with shape (n_fluid_cells,3), col headers = [rho, rho*v, e], where rho is fluid density, rho*v is momentum density, e is energy density. 
		
		t_start: double
			starting time of simulation
		
		t_end: double
			ending time of simulation

		x_start : double
			x position of first fluid cell

		x_end : double
			x position of last fluid cell
		
		n_cells: int 
			number of fluid cells excluding ghost cells
		
		boundary_pos: double
			x position of shock tube boundary at start_time
		
		boundary_left: array
			array of boundary conditions for left side of shocktube. Array expected to contain entries [rho, p, v], where rho is fluid density, p is pressure, v is velocity.
		
		boundary_right: array
			array of boundary conditions for right side of shocktube. Array expected to contain entries [rho, p, v], where rho is fluid density, p is pressure, v is velocity.
		
		boundary_condition: int
			int specifying how boundaries / ghost cells are updated. Choices = ['outflow','periodic', 'reflecting']
		
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

	// define independant and dependant variables along with step size vars
	double dx,dt,t,x; // grid_cell_size, time step, time and position vars;
	double v,rho,p,e; // velocity, density, pressure, energy density 
	dx = (x_end-x_start)/n_cells; 
	t = t_start;
	
	// define (N,3) array to keep track of cell fluxes at each time
	double (*F)[3] = malloc(sizeof(double[N+2][3]));

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
	for (int i=1; i<=n_cells; i++)
	{
		x = i*dx; 
		if (x<boundary_pos)
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
		// [STEP 1] set ghost / boundary cells
		set_ghosts(Q,boundary_condition);
		
		// [STEP 2] determine time step dt
		dt = get_cfl_dt(Q, dx, n_cells, gamma, 0.3);
		printf("%lf\n",dt);
		
		// [STEP 3] set fluxes within cells
		for (int i=0; i<n_cells+2; i++)
		{	
			p  = (gamma - 1) * Q[i][2];
			F[i][0] = Q[i][1];
			F[i][1] = pow(Q[i][1],2)/Q[i][0] + p;
			F[i][2] = Q[i][1] * (Q[i][2] + p_left) / Q[i][0];

			// printf("cell %i: \ninitial state:(%lf,%lf,%lf)\ninitial fluxes:(%lf,%lf,%lf)\n\n",i,Q[i][0],Q[i][1],Q[i][2],F[i][0],F[i][1],F[i][2]);
		}

		// FLUXES AND Q VECTORS ARE SET
		// [STEP 4] update cell state array Q using stepper
		//printf("t: %lf Q STATE 10 BEFORE update: (%lf,%lf,%lf), dx,dt = (%lf,%lf))\n",t,Q[10][0],Q[10][1],Q[10][2],dx,dt);
		stepper(Q,F,n_cells,dx,dt);
		//printf("Q STATE 10 AFTER update: (%lf,%lf,%lf)\n\n",Q[10][0],Q[10][1],Q[10][2]);
		// [STEP 5] update time
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
			fprintf(pf,"%lf,%lf,%lf,%lf,%lf\n", x, Q[i][0], (gamma - 1)*Q[i][2], Q[i][1]/Q[i][0], Q[i][2]);
		}
	}
}
