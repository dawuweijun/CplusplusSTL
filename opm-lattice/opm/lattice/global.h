
#ifndef HEADER
#define HEADER


/***********************************************************************/
/* Headerfile for LB code                                              */
/* D3Q19                                                               */
/***********************************************************************/

#define TOL 1.0e-12
#define PI 3.14159265359
#define SQRT2 1.41421356237
#define b 19
#define w0 1.0/3.0
#define w1 1.0/36.0
#define w2 1.0/18.0
#define c0 0.0 
#define c1 1.41421356
#define c2 1.0


/***********************************************************************/
/* GLOBAL PARAMETERS                                                   */
/***********************************************************************/

/* Effective dimension of latttice */
 const int Dim = 4; 

/* Total timesteps */
 const int T = 100000;  

/* Total density at each cite */
 const double RHO = 1.0; 
/* Total density at each cite */
// const double RHO1 = 0.8; 


/* Number of propagating sites */ 
 const int bm = 18; 

/* Number of rest sites */
 const int bc = 1;  

/* Rescale parameters time*/
 const double tau = 1.0;

/* Rescale parameters length*/
 const double h = 1.0; 

/* Default system size */
int Lx; 
int Ly;
int Lz; 

/* Default directional velocities */
 const double c_x[b] = {c2,-c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c2,-c2, c2,-c2, c0, c0, c0, c0, c0};
 const double c_y[b] = {c0, c0, c2,-c2, c0, c0, c2,-c2,-c2, c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c0}; 
 const double c_z[b] = {c0, c0, c0, c0, c2,-c2, c0, c0, c0, c0, c2,-c2,-c2, c2, c2,-c2,-c2, c2, c0}; 


/* Weight array */
 const double w[b] = {w2, w2, w2, w2, w2, w2,                          // Sublattice 1
			    w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1,  // Sublattice 2
			    w0};                                             // Sublattice 0 (rest particle)


/* Lattice viscosity */
 const double nuRed = 0.1667;  //red 
 const double nuBlue = 0.1667; //blue
 const double velMax = 0.1;  //max sound-speed


/* Global variables */
double ExForceConst;
double FluxForce;
double lambdaRed;
double lambdaBlue;


double ****f, ****g;      
double *N;

int ***Solid;

int ite;
int startTime;
int solidArea;
int nonSolidArea;

double J, Jfix;
double Jsum;


/*******************************************************************/
/* FUNCTIONS                                                       */
/*******************************************************************/

// Initializes phase distribution
// Input parameters: random seed, Size, phase arrays
// Manipulates Red and Blue arrays
void initializeFluid(int, int, int, double ****, double ****);

// Calculates the equlibrium distribution of phases for initialization
// Input: Density constant, position, local momentum
// Output: Local density
double feq(double, int, int, int, int, double, double, double);

// Runs simulations
// Input: Red and Blue arrays
void takeStep(double ****, double ****);

// Updates local populations
// Input: size, Red and Blue arrays
void update(int, int, int, double ****, double ****);

// Local collision step. Includes bounce-back, gradient calculation
// and phase separation.
// Input: Size, Red and Blue arrays
//void collisionStep(int, int, int, double ****, double ****);
void collisionStepSCRed(int, int, int, double ****, double ****);
void collisionStepSCBlue(int, int, int, double ****, double ****);


// Calculates local colour gradient for use in phase separation
// Input: Position and array to be manipulated
//void fcalc(int, int, int, double ****);
void fcalcSC(int, int, int,double ****, double *, double *, double ****);

// Local pseudo equilibrium distribution
// Input: Position, densities and momentums
// Output: pseudo equilibrium population
//double NiP(int, int, int, int, double, double, double, double, double);

// Same as above with zero velocity
//double NiP0(int, int, int, int, double, double);

double NiPSC(int, int, int, int, int, double, double, double, double, double);
//double NiPSCBlue(int, int, int, int, double, double, double, double, double);

// Calculation of surface tension
// Input: Position, total relaxation parameter, density and local gradient
// Output: Surface tension force
//double colorSurTen(int, int, int, int, double, double, double ****);

// Old phase separation function
//void colorRedistrib2(int, int, int, double,
//		     double *, double ****, double ****, double ****);

// Updated phase separation algorithm
//void colorRedistrib3(int, int, int, double, double,
//		     double *, double ****, double ****, double ****);

// Functions that take care of the boundary conditions of 
// of the streaming step
//void streamingStep(int, int, int, double ****, double ****);
void streamingStepBySwapping(int, int, int, double ****, double ****);

// Functions that write to file
void massMomentumCalc(double ****, double ****, int ***,char *, char *, char *, char *);
void pressureCalc(double ****, double ****, int ***);

// Allocation and de-allocation of memory
double ****create4DDouble(int, int, int, int);
int ***create3DInt(int, int, int);
void delete3DInt(int ***, int, int);
void delete4DDouble(double ****, int, int, int);

// Swapping to variables
void swap(double *, double *);

// Propagation functions. Propagates populations in velocity directions.
void propagationBySwapping(int, int, int, double ****);

/*******************************************************************/

#endif
