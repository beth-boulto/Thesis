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
#include <vector>
#include <cmath>
{
	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
Tdouble m=state[0];
Tdouble v=state[1];
Tdouble crit=state[2];
Tdouble u=control[0];
double phi=constants[0];
double muM=constants[1];
double alpha=constants[2];
double p=constants[3];
double r=constants[4];
double MT=constants[5];
double VT=constants[6];
double alp1=constants[7];
double alp2=constants[8];
double h12; 
double h13;
h12=r*p*pow((1-p),r)/(pow((1-p),(r+1)));
h13= r*(r+1)*p*p*pow((1-p),r)/(pow((1-p),r+2));

	state_dynamics[0] = phi*h12-alpha*v-(muM+u)*m;
state_dynamics[1]=phi*(h13+h12)+(muM+u)*m-2*(muM+u)*v;
state_dynamics[2]=u*u+alp1*(m-MT)*(m-MT)+alp2*(v-VT)*(v-VT);

}


