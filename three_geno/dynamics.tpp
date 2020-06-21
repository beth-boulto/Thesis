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
#include <cmath>
{
Tdouble u = control[0];
Tdouble ss=state[0];
Tdouble rr=state[2];
Tdouble sr=state[1];
Tdouble crit=state[3];
double beta=constants[0];
double muM=constants[1];
double DM=constants[2];
double psi1=constants[3];
double psi2=constants[4];
double psi3=constants[5];
double MT=constants[6];
double alp1=constants[7];
double PT=constants[8];
double alpr=constants[9];
double gamma=constants[10];

	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
	state_dynamics[0] = beta*(ss*ss+ss*sr+0.25*sr*sr)/(ss+sr+rr)-(muM+u)*ss-DM*(ss*ss+(1-psi1)*ss*sr+(1-psi2)*ss*rr);
state_dynamics[1] = beta*(2*ss*rr+ss*sr+0.5*sr*sr+rr*sr)/(ss+sr+rr)-(muM+gamma*u)*sr-DM*(sr*sr+(psi1)*ss*sr+(1-psi3)*sr*rr);
state_dynamics[2] = beta*(rr*rr+rr*sr+0.25*sr*sr)/(ss+sr+rr)-(muM)*rr-DM*(rr*rr+(psi3)*rr*sr+(psi2)*ss*rr);
state_dynamics[3]=u*u+alp1*(ss+rr+sr-MT)*(ss+rr+sr-MT)+alpr*((2*rr+sr)/(2*(rr+sr+ss))-PT)*((2*rr+sr)/(2*(rr+sr+ss))-PT);



}


