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
Tdouble Ms=state[0];
Tdouble Mr=state[1];
Tdouble Vs=state[2];
Tdouble Vr=state[3];
Tdouble crit=state[4];
double beta=constants[0];
double muM=constants[1];
double DM=constants[2];
double psi=constants[3];
double MT=constants[4];
double alp1=constants[6];
double VT=constants[5];
double alp2=constants[7];
double PT = constants[8];
double alpr=constants[9];
	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
	state_dynamics[0] = beta*Ms-(muM+u)*Ms-DM*(Ms*Ms+Vs+(1-psi)*Ms*Mr);
state_dynamics[1] = beta*Mr-(muM)*Mr-DM*(Mr*Mr+Vr+(psi)*Mr*Ms);
state_dynamics[2]=beta*Ms-(muM+u)*(2*Vs-Ms)-DM*(4*Vs*Vs/Ms+Vs*Ms-3*Vs-Ms*Ms+(1-psi)*(2*Mr*Vs-Mr*Ms));

state_dynamics[3]=beta*Mr-(muM)*(2*Vr-Mr)-DM*(4*Vr*Vr/Mr+Vr*Mr-3*Vr-Mr*Mr+(psi)*(2*Ms*Vr-Mr*Ms));
state_dynamics[4]=u*u+alp1*(Ms+Mr-MT)*(Ms+Mr-MT)+alp2*(Vs+Vr-VT)*(Vs+Vr-VT)+alpr*(Mr/(Ms+Mr)-PT)*(Mr/(Ms+Mr)-PT);


}


