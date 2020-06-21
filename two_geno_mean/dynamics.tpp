// Function for the dynamics of the problem
// dy/dt = dynamics(y,u,z,p)

// The following are the input and output available variables 
// for the dynamics of your optimal control problem.

// Input :
// time : current time (t)
// normalized_time: t renormalized in [0,1]
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization parameters
// constants : vector of constants

// Output :
// state_dynamics : vector giving the expression of the dynamic of each state variable.

// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])

// Tdouble variables correspond to values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during optimization use standard types (double, int, ...).

#include "header_dynamics"
{
Tdouble u = control[0];
Tdouble s=state[0];
Tdouble r=state[1];

Tdouble crit=state[2];
double beta=constants[0];
double muM=constants[1];
double DM=constants[2];
double psi=constants[3];
double MT=constants[4];
double alp1=constants[5];
double PT = constants[6];
double alpr = constants[7];
	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
	state_dynamics[0] = beta*s-(muM+u)*s-DM*(s*s+(1-psi)*s*r);
state_dynamics[1]=beta*r-(muM)*r-DM*(r*r+(psi)*s*r);
state_dynamics[2]=u*u+alp1*(s+r-MT)*(s+r-MT)+alpr*(r/(r+s)-PT)*(r/(r+s)-PT);
}


