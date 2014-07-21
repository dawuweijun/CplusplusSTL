# include <iostream>
# include <fstream>
# include <sstream>
# include <string>
# include <cstdio>
# include <cmath>
# include <cstdlib>
# include <sys/stat.h>
# include <sys/types.h>
# include "global.h"

/****************************************************************/
/* D3Q19 two-phase                                              */
/****************************************************************/

/****************************************************************/
/* MAIN METHOD                                                  */
/****************************************************************/

int main(int argc, char *argv[]) {

  /* Initialize *******************************/
//  long seed = 6781;
  int i;
  double nuRedScale;
  double nuBlueScale;

  /* Phase population distributions */
  double ****Red, ****Blue;

  Lx = 80; 
  Ly = 20;
  Lz = 20; 

  //Adds layers to ensure periodic BCs
 // Lx += (2 * layers);


  /* Allocate memory for global arrays. */
  Solid = create3DInt(Lz, Ly, Lx);
  Red = create4DDouble(Lz, Ly, Lx, b);
  Blue = create4DDouble(Lz, Ly, Lx, b);
  f = create4DDouble(Lz, Ly, Lx, b);
  N = new double[b];

  /**********************************************************/
  /* Initialization of solid matrix.                        */
  /*                                                        */
  /* Regular geometry for testing purposes                  */
  /*                                                        */
  /**********************************************************/

  

  // Scaled global constants 
  // If tau = 1 no scaling
  ExForceConst = tau*(double)Dim/(bm*c1*c1);
  nuRedScale = tau/(h*h)*nuRed;
  nuBlueScale = tau/(h*h)*nuBlue;
  
  // Relaxation parameters
  lambdaRed = (-2.0/(6*nuRedScale+1));
  lambdaBlue = (-2.0/(6*nuBlueScale+1));
  
  // Initial body-force
  FluxForce=0.0;

  // Initial volumetric flux
  J=0.0;
  Jfix = 0.0;

  startTime = 0;
  
  // timesteps
  ite = startTime;

  // Sets the initial phase populations
  initializeFluid(Lz,Ly,Lx,Red,Blue);

  // Calculates certain quantities and output to file
  // Not a part of the simulation
  massMomentumCalc(Red,Blue,Solid,"outputh","outputl","LBMflow","LBMflow");

  // Starts simulation
  takeStep(Red,Blue);
  
  // Cleaning up
  delete4DDouble(Red, Lz, Ly, Lx);
  delete4DDouble(Blue, Lz, Ly, Lx);
  delete4DDouble(f, Lz, Ly, Lx);
  delete3DInt(Solid, Lz, Ly);
  free(N);

  return 0;
}
/*******************************************************************/


/****************************************************************/
/* BASIC METHODES                                               */
/****************************************************************/

/****************************************************************/

int ***create3DInt(int nx, int ny, int nz)
{
    int ***p;   //declaration of p as: pointer-to-pointer-to-pointer of int 
    int x, y, z;
	
    p = new int**[nx]; //Allocate pointers for the nx 
    if (p != NULL)
    {
        for (x = 0; x < nx; x++)
        {
            p[x] = new int*[ny]; //Allocate pointers for the ny 
            if (p[x] == NULL)
            {
                printf("Memory allocation failed. Exiting....");
                //return (1);
            }
            else
            {
                for (y = 0; y < ny; y++)
                {
                    p[x][y] = new int[nz]; //Allocate pointers for the nz
                    if (p[x][y] == NULL)
                    {
                        printf("Memory allocation failed. Exiting....");
                        //return (1);
                    }
                }
            }
        }
    }
    else
    {
        printf("Memory allocation failed. Exiting....");
        //return (1);
    }
    return p;
}

/****************************************************************/
void delete3DInt(int ***p, int nx, int ny) 
{
  int x, y;

  	for (x = 0; x < nx; x++) 
	{
    	for(y = 0; y < ny; y++) 
		{
			delete[] p[nx][ny];
			delete[] p[nx];
	    }
  	}
  	delete[] p;
}

/****************************************************************/
double ****create4DDouble(int nx, int ny, int nz, int nu)
{
    double ****p; /*declaration of p as: pointer-to-pointer-to-pointer of int */
    int x, y, z;

    p = new double***[nx]; /*Allocate pointers for the nx */
    if (p != NULL)
    {
        for (x = 0; x < nx; x++)
        {
            p[x] = new double**[ny];/*Allocate pointers for the ny */
            if (p[x] == NULL)
            {
                printf("Memory allocation failed. Exiting....");
                //return (1);
            }
            else
            {
                for (y = 0; y < ny; y++)
                {
                    p[x][y] = new double*[nz]; /*Allocate pointers for the nz */
                    if (p[x][y] == NULL)
                    {
                        printf("Memory allocation failed. Exiting....");
                        //return (1);
                    }
                    else
                    {
                        for (z = 0; z < nz; z++)
                        {
                            p[x][y][z] = new double[nu]; /*Allocate pointers for the nz */
                            if (p[x][y][z] == NULL)
                            {
                                printf("Memory allocation failed. Exiting....");
                                //return (1);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        printf("Memory allocation failed. Exiting....");
        //return (1);
    }
    return p;

}

/****************************************************************/
void delete4DDouble(double ****p, int nx, int ny, int nz) 
{
  	int x, y, z;
  
  	for (x = 0; x < nx; x++) 
	{
    	for(y = 0; y < ny; y++) 
		{
      		for(z = 0; z < nz; z++) 
			{
				delete[] p[nx][ny][nz];
				delete[] p[nx][ny];
				delete[] p[nx];
      		}
    	}
  	}
  	delete[] p;
} 

/****************************************************************/

void initializeFluid(int Lz, int Ly, int Lx, 
		     double ****R, double ****B){
  
  int x, y, z, k;
  double ux,uy,uz;

 for(z = 0; z < Lz; z++) {
    for(y = 0; y < Ly; y++) {
      for(x = 0; x < Lx; x++) {
	Solid[z][y][x] = 1;
	

	  if(y < Ly - 1 && y > 0 && z < Lz - 1 && z > 0) {
	    Solid[z][y][x] = 0;
	  }

	
      }
    }
  }
  // According to bubble in the middle of a tube
  for(z = 0; z < Lz; z++) {
    for(y = 0; y < Ly; y++) {
      for(x = 0; x < Lx; x++) {
	if(Solid[z][y][x] != b) {
	  if(x < Lx/4 || x > 3*Lx/4) {
	    ux=0.0;
	    uy=0.0;
	    uz=0.0;
	    
	    for(k = 0; k < b; k++) {
	      
	      B[z][y][x][k] = feq(RHO,z,y,x,k,ux,uy,uz);
	      R[z][y][x][k] = 0.0;
	      
	    }
	  }
	  
	  else{
	    ux=0.0;
	    uy=0.0;
	    uz=0.0;
	    
	    for(k = 0; k < b; k++) {
	      R[z][y][x][k]= feq(RHO,z,y,x,k,ux,uy,uz);
	      B[z][y][x][k] = 0.0;
	      
	      
	    }
	  }
	}
      }
    }
  }
  

}//end initializeFluid()_________________________________________________

double feq(double density, int z, int i, 
		  int j, int k, double ux, double uy, double uz){
  double distribution, uc, u2;
  
  uc = ux*c_x[k] + uy*c_y[k] + uz*c_z[k];
  u2 = ux*ux + uy*uy + uz*uz;
  
  distribution = density*w[k]*(1.0 + 3*uc + 4.5*uc*uc - 1.5*u2);

  return distribution;

}//end feq()__________________________________________________

void massMomentumCalc(double ****R, double ****B, int ***Solid, char *directory1, char *directory2, char *filename1, char *filename2){
  
  int z, y, x, k, dir1, dir2;
  double mass = 0.0, redmass = 0.0, bluemass = 0.0, maxabsvelocity = -10.0;
  double uxtemp, uytemp, uztemp, densitytemp, temp1, temp2, val;
  int xmaxvel = -10, ymaxvel = -10, zmaxvel = -10;
  char dataFilename[255];
  FILE *dataFile1, *dataFile2;
  
  for(z = 0; z < Lz; z++) {
    for(y = 0; y < Ly; y++) {
      for(x = 0; x < Lx; x++) {
	uxtemp=0.0;
	uytemp=0.0;
	uztemp=0.0;
	densitytemp=0.0;
	
	for(k = 0; k < b; k++) {
	  mass += R[z][y][x][k] + B[z][y][x][k];
	  redmass += R[z][y][x][k];	  
	  bluemass += B[z][y][x][k];
	}

	if(Solid[z][y][x] == 0) {
	  for(k = 0; k < b; k++) {

	    densitytemp += R[z][y][x][k] + B[z][y][x][k];
	    uxtemp += (R[z][y][x][k] + B[z][y][x][k])*c_x[k];
	    uytemp += (R[z][y][x][k] + B[z][y][x][k])*c_y[k];
	    uztemp += (R[z][y][x][k] + B[z][y][x][k])*c_z[k];
	    
	  }
	  
	  uxtemp/=densitytemp;
	  uytemp/=densitytemp;
	  uztemp/=densitytemp;
	
	  if(sqrt(uxtemp*uxtemp+uytemp*uytemp+uztemp*uztemp) > maxabsvelocity){
	    
	    maxabsvelocity = sqrt(uxtemp*uxtemp+uytemp*uytemp+uztemp*uztemp);
	    xmaxvel = x;
	    ymaxvel = y;
	    zmaxvel = z;
	  }
	}
      }	 
    } 
  }
  
  
  ////////////////////////////
  //Dumping graphics to file//
  ////////////////////////////
  if(ite%10  == 0) {

#ifdef WIN32
    dir1 = mkdir(directory1);
    dir2 = mkdir(directory2);
#else
  //  dir = mkdir(directory,0777);
    dir1 = mkdir(directory1, 0777);
    dir2 = mkdir(directory2, 0777);
#endif
   if (dir1==0) printf("Error: Can't create output directory!\n");
    sprintf(dataFilename,"%s/%s_%07d.vti",directory1,filename1,ite);
    dataFile1 = fopen(dataFilename,"w");
    fprintf(dataFile1, "<?xml version=\"1.0\"?>\n");
    fprintf(dataFile1, "<!-- LBMflow v1.0.1, www.lbmflow.com -->\n");
    fprintf(dataFile1, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(dataFile1, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",Lx-1,Ly-1,Lz-1);
    fprintf(dataFile1, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Lx-1,Ly-1,Lz-1);
    fprintf(dataFile1, "    <PointData Scalars=\"scalars\">\n");
    fprintf(dataFile1, "       <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    
    for(z = 0; z < Lz; z++) {
      for(y = 0; y < Ly; y++) {
	for(x = 0; x < Lx; x++) {
	  temp1 = 0.0;
	  for(k = 0; k < b; k++) {
	    temp1 += R[z][y][x][k];
	    
	  }

	  val = temp1 / (RHO + 1.0);
	  if(Solid[z][y][x] != 0) val = 0.;

	fprintf(dataFile1, "%f ", val);

	}
	fprintf(dataFile1, "\n");
      }
    }
fprintf(dataFile1, "      </DataArray>\n");

    //write pressure

    fprintf(dataFile1, "    </PointData>\n");

    fprintf(dataFile1, "    <CellData>\n");
    fprintf(dataFile1, "    </CellData>\n");
    fprintf(dataFile1, "  </Piece>\n");
    fprintf(dataFile1, "  </ImageData>\n");

    fprintf(dataFile1, "</VTKFile>\n");
    fclose(dataFile1);

if (dir2==0) printf("Error: Can't create output directory!\n");
    sprintf(dataFilename,"%s/%s_%07d.vti",directory2,filename2,ite);
    dataFile2 = fopen(dataFilename,"w");
    fprintf(dataFile2, "<?xml version=\"1.0\"?>\n");
    fprintf(dataFile2, "<!-- LBMflow v1.0.1, www.lbmflow.com -->\n");
    fprintf(dataFile2, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(dataFile2, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",Lx-1,Ly-1,Lz-1);
    fprintf(dataFile2, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Lx-1,Ly-1,Lz-1);
    fprintf(dataFile2, "    <PointData Scalars=\"scalars\">\n");
    fprintf(dataFile2, "       <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for(z = 0; z < Lz; z++) {
      for(y = 0; y < Ly; y++) {
	for(x = 0; x < Lx; x++) {
	  temp2 = 0.0;
	  for(k = 0; k < b; k++) {
	    temp2 += B[z][y][x][k];
	    
	  }
	  val = temp2 / (RHO + 1.0);
	  if(Solid[z][y][x] != 0) val = 0.;
	fprintf(dataFile2, "%f ", val);

	}
	fprintf(dataFile2, "\n");
      }
    }
fprintf(dataFile2, "      </DataArray>\n");

    //write pressure

    fprintf(dataFile2, "    </PointData>\n");

    fprintf(dataFile2, "    <CellData>\n");
    fprintf(dataFile2, "    </CellData>\n");
    fprintf(dataFile2, "  </Piece>\n");
    fprintf(dataFile2, "  </ImageData>\n");

    fprintf(dataFile2, "</VTKFile>\n");
    fclose(dataFile2);
  }




  printf("Total mass = %f\n", mass);  
  printf("Red mass = %f\n", redmass);
  printf("Blue mass = %f\n", bluemass);
  printf("F0 = %f\n", FluxForce);
  printf("Maximal abs. velocity = %f, (%f, %f %f)\n", maxabsvelocity,  
	 uxtemp, uztemp, uytemp); 
    
}/* end massMomentumCalc()____________________________________________ */
void takeStep(double ****R, double ****B) {

  while(ite < T) {

    // updates phase populations
    update(Lz, Ly, Lx, R, B);

    // Writes to file
    if(ite%100 == 0) {
      massMomentumCalc(R, B, Solid,"outputh", "outputl", "LBMflow", "LBMflow");
      pressureCalc(R, B, Solid);
    }
    
    ite++;  /* Updates total timesteps ellapsed ****************/
    
  }

}
void pressureCalc(double ****R, double ****B, int ***Solid) {
  
  int z, y, x, k, nonSolidPressureRed, nonSolidPressureBlue;
  double pressureRed = 0.0;
  double pressureBlue = 0.0;
  double Jred, Jblue, densityRed, densityBlue, density, densityRedTemp, densityBlueTemp;
  double uxtemp, uytemp, uztemp, densitytemp,  xpos, ux, force, momentum;
  FILE *outFile1;
  char fname[30];


  sprintf(fname, "pressureFile_%i.dat", ite);
  outFile1 = fopen(fname, "w");

  for(x = 0; x < Lx; x++) {
    
    if(x%2 == 0) {
      if(nonSolidPressureRed > 0) {
	pressureRed /= (double) nonSolidPressureRed;
      }
      
      if(nonSolidPressureBlue > 0) {
	pressureBlue /= (double) nonSolidPressureBlue;
      }
      fprintf(outFile1, "%i %e %e\n", x, pressureRed, pressureBlue);
      pressureRed = 0.0;
      pressureBlue = 0.0;
      nonSolidPressureRed = 0;
      nonSolidPressureBlue = 0;
      
    }

    for(y = 0; y < Ly; y++) {
      for(z = 0; z < Lz; z++) {
	uxtemp=0.0;
	uytemp=0.0;
	uztemp=0.0;
	densitytemp=0.0;
	densityRedTemp = 0.0;
	densityBlueTemp = 0.0;
	
	if(Solid[z][y][x] == 0) {
	  for(k = 0; k < b; k++) {
	    densitytemp += R[z][y][x][k] + B[z][y][x][k];
	    densityRedTemp += R[z][y][x][k];
	    densityBlueTemp += B[z][y][x][k];
	    uxtemp += (R[z][y][x][k] + B[z][y][x][k])*c_x[k];
	    uytemp += (R[z][y][x][k] + B[z][y][x][k])*c_y[k];
	    uztemp += (R[z][y][x][k] + B[z][y][x][k])*c_z[k];
	    force += fabs(ExForceConst * w[k] * FluxForce * RHO * c_x[k]);
	    
	  }
	  

	  Jred += uxtemp*densityRedTemp/densitytemp;
	  Jblue += uxtemp*densityBlueTemp/densitytemp;

	  density += densitytemp;
	  densityRed += densityRedTemp;
	  densityBlue += densityBlueTemp;

	  momentum += uxtemp;

	  uxtemp/=densitytemp;
	  uytemp/=densitytemp;
	  uztemp/=densitytemp;

	  ux += uxtemp;
	  
	  if(densityRedTemp > 0.8*densitytemp) {
	    nonSolidPressureRed++;
	    pressureRed += 0.33*densitytemp*(1 - (uxtemp*uxtemp + uytemp*uytemp + uztemp*uztemp))
	      + FluxForce*densitytemp*(Lx - x);
	  }
	  
	  if(densityBlueTemp > 0.8*densitytemp) {
	    nonSolidPressureBlue++;
	    pressureBlue += 0.33*densitytemp*(1 - (uxtemp*uxtemp + uytemp*uytemp + uztemp*uztemp))
	      + FluxForce*densitytemp*(Lx - x);
	  }
	}
	
      }
    }
  }
  

  fclose(outFile1);

  
}
	
void update(int Lz, int Ly, int Lx, double ****R, double ****B) {

  int z, y, x, k;

  // Collision
  // Chose if you want Color gradient or Shan-Chen

  // Color gradient
//  collisionStep(Lz,Ly,Lx,R,B);

  //Shan-Chen
   collisionStepSCRed(Lz,Ly,Lx,R,B);
   collisionStepSCBlue(Lz,Ly,Lx,R,B);
  
  // Streaming
  // Swapping trick saves you memory
  streamingStepBySwapping(Lz,Ly,Lx,R,B);

}
void collisionStepSCBlue(int lz, int ly, int lx, 
		   double ****R,
		   double ****B) {

  int i, z, y, x, k, nonSolid;
  double lambda, RedDensity, BlueDensity, ux, uy, uz, uxtilde, uytilde, uztilde;
  double Rtemp, Btemp;

  /* Shan-Chen "potential" strengths */
  /* Care must be taken when these parameters are set! */
  const double G0 = 0.0;
  const double G1 = -0.2;
  const double G2 = -0.1;
  const double GS = 0.01;

  double GFF[19] = {G1, G1, G1, G1, G1, G1, 
		    G2, G2, G2, G2, G2, G2,
		    G2, G2, G2, G2, G2, G2,
		    G0};
  
  double GFS[19] = {GS, GS, GS, GS, GS, GS, 
		    GS, GS, GS, GS, GS, GS,
		    GS, GS, GS, GS, GS, GS, 
		    G0};
  
  /****************************************************************/
   
  for(z = 0; z < lz; z++) {
    for(y = 0; y < ly; y++) {
      for(x = 0; x < lx; x++) {
	
	if(Solid[z][y][x] == 0) {
	  RedDensity = 0.0;
	  BlueDensity = 0.0;
	  for(k = 0; k < b; k++) {
	    RedDensity += R[z][y][x][k];
	    BlueDensity += B[z][y][x][k];
	    
	  }
	  
	  for(k = 0; k < b; k++) {
	    f[z][y][x][k] = RedDensity;
	  }
	  
	}
	
	else if(Solid[z][y][x] == 1) {
	  for(k = 0; k < b; k++) {
	    f[z][y][x][k] = -10.0;
	  }
	}
	
      }
    }
  }

  /* Propagates the densities to neighbors */  
  propagationBySwapping(lz, ly, lx, f);

  /* Calculates the SC potentials */
  fcalcSC(lz, ly, lx, f, GFF, GFS, f);

  /****************************************************************/  
  /* Bounce-back                                                  */
  /****************************************************************/  
    double RedAfter[b] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  double BlueAfter[b] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(z = 0; z < lz; z++) {
    for(y = 0; y < ly; y++) {
      for(x = 0; x < lx; x++) {
	
	if(Solid[z][y][x] == 1) {
	  
	  for(k = 0; k < bm; k++) {
	    if(k%2 == 0) {
	      BlueAfter[k] = B[z][y][x][k + 1]; 
	    }
	    
	    else {
	      BlueAfter[k] = B[z][y][x][k - 1]; 
	    }
	  }

	  for(k = 0; k < bm; k++) {	
	    B[z][y][x][k] = BlueAfter[k]; 
	  }
	  

	}
	
	else if(Solid[z][y][x] == 0){
	  
	  RedDensity=0.0; 
	  BlueDensity=0.0;
	  
	  ux = 0.0;
	  uy = 0.0;
	  uz = 0.0;
	  

	  for(k = 0; k < b; k++) {
	    Rtemp = R[z][y][x][k];
	    Btemp = B[z][y][x][k];
	    RedDensity += Rtemp;
	    BlueDensity += Btemp;
	    
	    ux += ( Rtemp + Btemp ) * c_x[k]
	      - (1.0/lambdaBlue) * (fabs(Btemp*FluxForce * c_x[k]));
	    
	    uy += ( Rtemp + Btemp ) * c_y[k];
	    
	    uz += ( Rtemp + Btemp ) * c_z[k];
	  }
	  

	  /* SC potential added to momentum */	  
	  ux += -(1.0/lambdaBlue) * BlueDensity * f[z][y][x][0];
	  uy += -(1.0/lambdaBlue) * BlueDensity * f[z][y][x][1];
	  uz += -(1.0/lambdaBlue) * BlueDensity * f[z][y][x][2];
	  
	  
	  
	  if((RedDensity + BlueDensity) <= 0.0)
	    lambda = 0.0;
	  else{
	    lambda = (RedDensity * lambdaRed + BlueDensity * lambdaBlue)
	      /(RedDensity + BlueDensity);
	  }
	  
	  for(k = 0; k < bm; k++) {
	    B[z][y][x][k] += lambdaBlue*(B[z][y][x][k] - NiPSC(0,z,y,x,k,RedDensity,BlueDensity,ux,uy,uz));
	  }
	  
	  B[z][y][x][b-1] += lambdaBlue*(B[z][y][x][b-1] - NiPSC(0,z,y,x,b-1,RedDensity,BlueDensity,ux,uy,uz));
	  

	}
      }
    }
  }
  
}
void collisionStepSCRed(int lz, int ly, int lx, 
		   double ****R,
		   double ****B) {

  int i, z, y, x, k, nonSolid;
  double lambda, RedDensity, BlueDensity, ux, uy, uz, uxtilde, uytilde, uztilde;
  double Rtemp, Btemp;

  /* Shan-Chen "potential" strengths */
  /* Care must be taken when these parameters are set! */
  const double G0 = 0.0;
  const double G1 = -0.2;
  const double G2 = -0.1;
  const double GS = -0.01;

  double GFF[19] = {G1, G1, G1, G1, G1, G1, 
		    G2, G2, G2, G2, G2, G2,
		    G2, G2, G2, G2, G2, G2,
		    G0};
  
  double GFS[19] = {GS, GS, GS, GS, GS, GS, 
		    GS, GS, GS, GS, GS, GS,
		    GS, GS, GS, GS, GS, GS, 
		    G0};
  
  /****************************************************************/
   
  for(z = 0; z < lz; z++) {
    for(y = 0; y < ly; y++) {
      for(x = 0; x < lx; x++) {
	
	if(Solid[z][y][x] == 0) {
	  RedDensity = 0.0;
	  BlueDensity = 0.0;
	  for(k = 0; k < b; k++) {
	    RedDensity += R[z][y][x][k];
	    BlueDensity += B[z][y][x][k];
	    
	  }
	  
	  for(k = 0; k < b; k++) {
	    f[z][y][x][k] = BlueDensity;
	  }
	  
	}
	
	else if(Solid[z][y][x] == 1) {
	  for(k = 0; k < b; k++) {
	    f[z][y][x][k] = -10.0;
	  }
	}
	
      }
    }
  }
  
  /* Propagates the densities to neighbors */
  propagationBySwapping(lz, ly, lx, f);

  /* Calculates the SC potentials */
  fcalcSC(lz, ly, lx, f, GFF, GFS, f);
  
  /****************************************************************/  
  /* Bounce-back                                                  */
  /****************************************************************/  
  
  double RedAfter[b] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  double BlueAfter[b] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(z = 0; z < lz; z++) {
    for(y = 0; y < ly; y++) {
      for(x = 0; x < lx; x++) {
	
	if(Solid[z][y][x] == 1) {
	  
	  for(k = 0; k < bm; k++) {
	    if(k%2 == 0) {
	      RedAfter[k] = R[z][y][x][k + 1];
	    }
	    
	    else {
	      RedAfter[k] = R[z][y][x][k - 1];
	    }
	  }

	  for(k = 0; k < bm; k++) {	
	    R[z][y][x][k] = RedAfter[k];
	  }
	  

	}
	
	else if(Solid[z][y][x] == 0){
	  
	  RedDensity=0.0; 
	  BlueDensity=0.0;
	  
	  ux = 0.0;
	  uy = 0.0;
	  uz = 0.0;

	  for(k = 0; k < b; k++) {
	    Rtemp = R[z][y][x][k];
	    Btemp = B[z][y][x][k];
	    RedDensity += Rtemp;
	    BlueDensity += Btemp;
	    
	    ux += ( Rtemp + Btemp ) * c_x[k]
	      - (1.0/lambdaRed) * (fabs(Rtemp*FluxForce*c_x[k]));
	    
	    uy += ( Rtemp + Btemp ) * c_y[k];
	    
	    uz += ( Rtemp + Btemp ) * c_z[k];
	  }
	  

	  /* SC potential added to momentum */
	  ux += -(1.0/lambdaRed) * RedDensity * f[z][y][x][0];
	  uy += -(1.0/lambdaRed) * RedDensity * f[z][y][x][1];
	  uz += -(1.0/lambdaRed) * RedDensity * f[z][y][x][2];
	  
	  
	  
	  if((RedDensity + BlueDensity) <= 0.0)
	    lambda = 0.0;
	  else{
	    lambda = (RedDensity * lambdaRed + BlueDensity * lambdaBlue)
	      /(RedDensity + BlueDensity);
	  }
	  
	  for(k = 0; k < bm; k++) {
	    R[z][y][x][k] += lambdaRed*(R[z][y][x][k] - NiPSC(1,z,y,x,k,RedDensity,BlueDensity,ux,uy,uz));
	  }
	  
	  R[z][y][x][b-1] += lambdaRed*(R[z][y][x][b-1] - NiPSC(1,z,y,x,b-1,RedDensity,BlueDensity,ux,uy,uz));
	  

	}
      }
    }
  }
  
}

void streamingStepBySwapping(int Lz, int Ly, int Lx, 
		   double ****R,
		   double ****B) {
  
  int z, y, x, k;
  double redDensityLx = 0.0, redBackDensityLx = 0.0;

  // Propagation step
   propagationBySwapping(Lz,Ly,Lx,R);
  propagationBySwapping(Lz,Ly,Lx,B);
}

void propagationBySwapping(int z, int y, int x, double ****Lattice) {

  double temp;
 int xtop, ytop, ztop, xbot, ybot, zbot;

  // First swap locally
  for(z = 0; z < Lz; z++) {

    for(y = 0; y < Ly; y++) {

      for(x = 0; x < Lx; x++) {

	if(Solid[z][y][x] != 2) {
	
	  swap(&Lattice[z][y][x][0], &Lattice[z][y][x][1]);
	  swap(&Lattice[z][y][x][2], &Lattice[z][y][x][3]);
	  swap(&Lattice[z][y][x][4], &Lattice[z][y][x][5]);
	  swap(&Lattice[z][y][x][6], &Lattice[z][y][x][7]);
	  swap(&Lattice[z][y][x][8], &Lattice[z][y][x][9]);
	  swap(&Lattice[z][y][x][10], &Lattice[z][y][x][11]);
	  swap(&Lattice[z][y][x][12], &Lattice[z][y][x][13]);
	  swap(&Lattice[z][y][x][14], &Lattice[z][y][x][15]);
	  swap(&Lattice[z][y][x][16], &Lattice[z][y][x][17]);
	}
      }
    }
  }	  

  // Swap between neigbours
  for(z = 0; z < Lz; z++) {
		ztop = (z >= Lz-1 ? z+1-Lz : z+1);
	zbot = (z <= 0 ? z-1+Lz : z-1);
    for(y = 0; y < Ly; y++) {
	ytop = (y >= Ly-1 ? y+1-Ly : y+1);
	ybot = (y <= 0 ? y-1+Ly : y-1);
      for(x = 0; x < Lx; x++) {
	xtop = (x >= Lx-1 ? x+1-Lx : x+1);
	xbot = (x <= 0 ? x-1+Lx : x-1);
	if(Solid[z][y][x] != 2) {

	  /***********************************************************/

 swap(&Lattice[z][y][x][0], &Lattice[z][y][xbot][1]);
	  swap(&Lattice[z][y][x][2], &Lattice[z][ybot][x][3]);
	  swap(&Lattice[z][y][x][4], &Lattice[zbot][y][x][5]);
	  swap(&Lattice[z][y][x][6], &Lattice[z][ybot][xbot][7]);
	  swap(&Lattice[z][y][x][8], &Lattice[z][ytop][xbot][9]);
	  swap(&Lattice[z][y][x][10], &Lattice[zbot][y][xbot][11]);
	  swap(&Lattice[z][y][x][12], &Lattice[ztop][y][xbot][13]);
	  swap(&Lattice[z][y][x][14], &Lattice[zbot][ybot][x][15]);
	  swap(&Lattice[z][y][x][16], &Lattice[ztop][ybot][x][17]);

	}
      }
    }
  }
  
}/* END FUNCTION */

void swap(double *val1, double *val2) {

  double temp;
  
  temp = *val1;
  *val1 = *val2;
  *val2 = temp;
    
}/* END FUNCTION */

double NiPSC(int fluid, int z, int y, int x, int k,
	   double Rdens, double Bdens, double ux, double uy, double uz) {
  int l;
  double Nip, uc, u2, rho = 0.0;

  rho = Rdens + Bdens;

  if(rho <= 0.0) {
    rho=0.0;
  }
  else{ 
    ux /= rho;
    uy /= rho;
    uz /= rho;
  }
  
  uc = ux*c_x[k] + uy*c_y[k] + uz*c_z[k];
  u2 = ux*ux + uy*uy + uz*uz;

if (fluid == 0)
{
  Nip = Bdens*w[k]*(1+3*uc+4.5*uc*uc-1.5*u2);
}
if (fluid == 1)
{
  Nip = Rdens*w[k]*(1+3*uc+4.5*uc*uc-1.5*u2);
}  
  return Nip;

}


void fcalcSC(int Lz, int Ly, int Lx, 
	     double ****f, double *gff, double *gfs, double ****F) {
  int z, y, x, l;
  double fx, fy, fz;
  
  for(z = 0; z < Lz; z++) {
    for(y = 0; y < Ly; y++) {
      for(x = 0; x < Lx; x++) {
	if(Solid[z][y][x] != 2) {
	  fx=0.0;
	  fy=0.0;
	  fz=0.0;

	  for(l = 0; l < b; l++) {

	    if(f[z][y][x][l] > 0.0) {
	      fx += gff[l]*c_x[l]*f[z][y][x][l];
	      fy += gff[l]*c_y[l]*f[z][y][x][l];
	      fz += gff[l]*c_z[l]*f[z][y][x][l];
	    }

	    else if(f[z][y][x][l] < -9.9) {
	      fx += gfs[l]*c_x[l];
	      fy += gfs[l]*c_y[l];
	      fz += gfs[l]*c_z[l];
	    }
	  }	
	  
	  F[z][y][x][0] = -fx;
	  F[z][y][x][1] = -fy;
	  F[z][y][x][2] = -fz;
	}
      }
    }
  }
}



