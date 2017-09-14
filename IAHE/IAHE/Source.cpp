#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <climits>
#include <iomanip>
#include <windows.h>
#include <sstream>
#include <cstring>
#include "stopwatch.h"
using namespace std;


int** allbest;
int** allswarm;
int** allbaseswarm;
long*** deltaArray;
int* allparticalparentindex;
double*** allvelocity;
long* allfittness;
long* allbestfittness;
float inertia;
float cof1;
float cof2;
float cof3;
float maxvelocity;
float minvelocity;
int height;
int degree;
long systembestfitness;
long* systembestfitnesslist;
int* iterationnumber;
int** systembestresultlist;
int* systembestresult;
int particlesize;
int* systemfunctioncounter;
int* systemParticleIndex;
double* PEAValues;
int* GACount;
int* PACount;
int functioncounter;
bool geneticalgorithm;
bool hybridgeneticalgorithm;
bool adaptivehybridgeneticalgorithm;
bool HPSO;
long bestfitnessknown;
int RTSfailures;
bool bestfittnessupdated;
int** B;
int** A;
int variablesize;
long** heuristic;
int systemIterator;
long int seed;
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
long maxiterations;
typedef int*   type_vector;
typedef long** type_matrix;
const long infinite = 2147483647;
Stopwatch swatch;
int step;
int mean;
int stdv;
int upperlimittenure;
int lowerlimittenure;
int upperlimitrestart;
int lowerlimitrestart;
int tenure_min;
int tenure_max;
int iteration;
int** CurrentGeneration;
int** NewPopulation;
float MutationProb;
float CrossoverProb;
long* CurrentGenFitness;
long* NewPopFitness;
int** SelectedParents;
int** Children;
int FPR;
int FNR;
double* systemtimeinSec;
int failurelimit;
int GAfailurelimit;
int* ParentIndexes;
long* SelectedFitness;
int* SelectedIndexes;
long FitnessComp[4];
int Index[4];
double probabilityGA;
bool particalechosenGA;
int GAAlgorithm;
int HPSOAlgorithm;
int SizeData;

/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
long compute_delta(int index, int i, int j)
{
	long d; int k;
	d = (A[i][i] - A[j][j])*(B[allbaseswarm[index][j]][allbaseswarm[index][j]] - B[allbaseswarm[index][i]][allbaseswarm[index][i]]) +
		(A[i][j] - A[j][i])*(B[allbaseswarm[index][j]][allbaseswarm[index][i]] - B[allbaseswarm[index][i]][allbaseswarm[index][j]]);
	for (k = 0; k < variablesize; k = k + 1) if (k != i && k != j)
		d = d + (A[k][i] - A[k][j])*(B[allbaseswarm[index][k]][allbaseswarm[index][j]] - B[allbaseswarm[index][k]][allbaseswarm[index][i]]) +
		(A[i][k] - A[j][k])*(B[allbaseswarm[index][j]][allbaseswarm[index][k]] - B[allbaseswarm[index][i]][allbaseswarm[index][k]]);
	return(d);
}

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
long compute_delta_part(int index,
	int i, int j, int r, int s)
{
	return(deltaArray[index][i][j] + (A[r][i] - A[r][j] + A[s][j] - A[s][i])*
		(B[allbaseswarm[index][s]][allbaseswarm[index][i]] - B[allbaseswarm[index][s]][allbaseswarm[index][j]] + B[allbaseswarm[index][r]][allbaseswarm[index][j]] - B[allbaseswarm[index][r]][allbaseswarm[index][i]]) +
		(A[i][r] - A[j][r] + A[j][s] - A[i][s])*
		(B[allbaseswarm[index][i]][allbaseswarm[index][s]] - B[allbaseswarm[index][j]][allbaseswarm[index][s]] + B[allbaseswarm[index][j]][allbaseswarm[index][r]] - B[allbaseswarm[index][i]][allbaseswarm[index][r]]));
}

bool ifContains(int* ind, int bit, int IndSize)
{
	for(int i = 0; i < IndSize; i ++)
		if(ind[i] == bit)
			return true;
	return false;
}

void Dispose()
{
	delete[] allbestfittness;
	delete[] allfittness;
	delete[] systembestfitnesslist;
	delete[] allparticalparentindex;
	delete[] iterationnumber;
	delete[] systembestresult;
	delete[] systemfunctioncounter;
	delete[] systemParticleIndex;
	delete[] systemtimeinSec;
	
	for (int i=0; i < particlesize; i = i+1) 
		delete[] allbest[i];
	delete[] allbest;

	for (int i=0; i < particlesize; i = i+1) 
		delete[] allbaseswarm[i];
	delete[] allbaseswarm;

	for (int i=0; i < particlesize; i = i+1) 
		delete[] allswarm[i];
	delete[] allswarm;

	for (int i=0; i < particlesize; i = i+1) 
	{
		for( int j=0; j < variablesize;j++)
			delete[] allvelocity[i][j];
		delete[] allvelocity[i];
	}
	delete[] allvelocity;

	for (int i=0; i < particlesize; i = i+1) 
	{
		for( int j=0; j < variablesize;j++)
			delete[] deltaArray[i][j];
		delete[] deltaArray[i];
	}
	delete[] deltaArray;

	for (int i=0; i < SizeData; i = i+1) 
		delete[] systembestresultlist[i];
	delete[] systembestresultlist;

	for (int i=0; i < variablesize; i = i+1) 
		delete[] B[i];
	delete[] B;

	for (int i=0; i < variablesize; i = i+1) 
		delete[] A[i];
	delete[] A;
	
	for (int i=0; i < variablesize; i = i+1) 
		delete[] heuristic[i];
	delete[] heuristic;
}


//random number between 1
double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

//Gaussian Generator
double generate_gaussian (double mean, double std)
{
        
        double gw, gy1 ; 
        double gx1, gx2 ; 
        double gu ;       
      
         do {
                gx1 = 2.0 * ran01 (&seed) - 1.0 ; 
                gx2 = 2.0 * ran01 (&seed) - 1.0 ; 
                gw = gx1 * gx1 + gx2 * gx2 ; 
        } while (gw >= 1.0) ; 
                      
        gw = sqrt ( (-2 * log (gw) ) / gw ) ;
        gy1 = (gx1 * gw) * std + mean  ;
          
        return gy1 ; 
}

/*********** return an integer between low and high *************/
int unif(int low, int high)
{return low + (int)((double)(high - low + 1) * ran01(&seed)) ;}


void intializeSwarm()
{
	
	height = 2;
	degree = 3;
	inertia = 0.729;
	cof1 = 1.494;
	cof2 = 1.494;
	maxvelocity = 2;
	minvelocity = -2;
	particlesize =0;
	for( int i=0;i<=height;i++)
	{
		particlesize += pow(degree,i);
	}
	mean = 0;
	stdv = 1;
	allbest = new int *[particlesize];
	allfittness = new long [particlesize];
	allbestfittness = new long [particlesize];
	allswarm = new int *[particlesize];
	allbaseswarm = new int *[particlesize];
	for(int i=0;i<particlesize;i++)
	{
		allswarm[i]  = new int [variablesize];
		allbest[i]  = new int [variablesize];
		allbaseswarm[i]  = new int[variablesize];
	}
	allvelocity = new double **[particlesize];
	deltaArray = new long **[particlesize];
	for(int i=0;i<particlesize;i++)
	{
		allvelocity[i] = new double *[variablesize];
		deltaArray[i] = new long *[variablesize];
		for(int j=0;j< variablesize;j++)
		{
			allvelocity[i][j] = new double [variablesize];
			deltaArray[i][j] = new long[variablesize];
		}
	}
	allparticalparentindex = new int [particlesize];
	functioncounter = 0;
	systemfunctioncounter = new int[SizeData];
	systembestresult = new int [SizeData];
	systemParticleIndex = new int[SizeData];
	systembestresultlist = new int *[SizeData];
	systemtimeinSec = new double[SizeData];
	for( int i=0; i < SizeData ;i++)
	{
		systembestresultlist[i] = new int [ variablesize];
	}
	systembestfitnesslist = new long [SizeData];
	systemIterator = 0;
}

//read the problem files and add them to A,B matrix and creates heuristic matrix 
void readFile(string filename)
{
	// initialize stream reader
	ifstream read;
	read.open(filename.c_str(),ios::in);
	
	//read the problem size
	read >> variablesize;
	//initialize the pointers used;
	A = new int *[variablesize];
	B = new int *[variablesize];
	int* AArray = new int [variablesize];
	int* BArray = new int [variablesize];
	//initialize the 3  2d arrays 
	for(int i=0 ; i<variablesize ; i++)
	{
		A[i] = new int [variablesize];
		B[i] = new int [variablesize];
	}
	
	for (int i = 0; i < variablesize; i = i+1) 
		for (int j = 0; j < variablesize; j = j+1)
			read >> A[i][j];

    for (int i = 0; i < variablesize; i ++) 
		for (int j = 0; j < variablesize; j ++)
			read >> B[i][j];
	
	
	delete[] AArray;
	delete[] BArray;
	 read.close();
}

void transpose(int &a, int &b) 
{
	long temp = a; 
	a = b; 
	b = temp;
}

//generate random solution
void gererateRandom(int index)
{
	for (int i = 0; i < variablesize; i++)
    {
            allswarm[index][i] = i;
	}
     for (int i = 0; i < variablesize; i++) 
     {
         transpose(allswarm[index][i], allswarm[index][unif(i, variablesize - 1)]);
     }
}

void generateRandomVelocity(int index)
 {
	 for(int i=0;i<variablesize;i++)
	 {
		 for(int j=0;j<variablesize;j++)
	     {
			 allvelocity[index][i][j] = ran01(&seed);
		 }
	 }
 }

int evaluate(int index)
{
	int total = 0;
	for (int i = 0; i < variablesize; i++)
    {
        for (int j = 0; j < variablesize; j++)
         {
			 total += (A[i][j]) * (B[allswarm[index][i]][allswarm[index][j]]); 
         }
	}
    return total;
}

void intializeheirarchy()
{
      int bestpartical = -1;
       int parentCounter = 0;
	   for (int i = 0; i < particlesize; i++)
       {

           if (i == 0)
           {
			   allparticalparentindex[i] = bestpartical;
            }
            else
            {
				for (int l = 0; l < degree; l++)
                {
                        allparticalparentindex[i] = parentCounter;
                        i++;
                }
                i--;
                parentCounter++;
             }

       }
}
 
void calculateFittness(int index,bool type)
 {
	 functioncounter++;
	 long fitnessparticle = allfittness[index];
	 for(int i_retained=0;i_retained<variablesize;i_retained++)
	 {
		 if(allbaseswarm[index][i_retained] != allswarm[index][i_retained])
		 {
			 for(int j_retained = i_retained + 1; j_retained<variablesize; j_retained++)
			 {
				 if(allbaseswarm[index][j_retained] == allswarm[index][i_retained])
				 {
					 transpose(allbaseswarm[index][i_retained],allbaseswarm[index][j_retained]);
					 fitnessparticle += deltaArray[index][i_retained][j_retained];
					 for (int i = 0; i < variablesize -1; i = i+1) 
					 {	
						for (int j = i+1; j < variablesize; j = j+1)
						{	
							if (i != i_retained && i != j_retained && 
								j != i_retained && j != j_retained)
								 {
									 deltaArray[index][i][j] = compute_delta_part(index,i, j, i_retained, j_retained);
								}
								else
								{
									 deltaArray[index][i][j] = compute_delta(index, i , j);
								}
						}
					 }
					 break;
				 }
			 }
		 }
	 }
	 if(fitnessparticle < allbestfittness[index])
	 {
		 for(int j = 0 ; j< variablesize ; j++)
		 {
			allbest[index][j] = allswarm[index][j];
		 }
		 if(type)
			 GAAlgorithm++;
		 else
			 HPSOAlgorithm++;
		 allbestfittness[index] = fitnessparticle;
		 if(fitnessparticle < systembestfitness)
		 {
			 bestfittnessupdated = true;
			 systembestfitness = fitnessparticle;
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresult[j] = allswarm[index][j];
			 }
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresultlist[systemIterator][j] = allswarm[index][j];
			 }
			 systemfunctioncounter[systemIterator] = iteration;
			 systemParticleIndex[systemIterator] = index;
			 systembestfitnesslist[systemIterator] = fitnessparticle;
			 systemtimeinSec[systemIterator] = swatch.get_time_so_far("FullRun");
			 systemIterator ++; 
		 }
	 }
	 allfittness[index] = fitnessparticle;
 }


/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
long compute_delta(type_vector & p, int i, int j)
{
	long d; int k;
	d = (A[i][i] - A[j][j])*(B[p[j]][p[j]] - B[p[i]][p[i]]) +
		(A[i][j] - A[j][i])*(B[p[j]][p[i]] - B[p[i]][p[j]]);
	for (k = 0; k < variablesize; k = k + 1) if (k != i && k != j)
		d = d + (A[k][i] - A[k][j])*(B[p[k]][p[j]] - B[p[k]][p[i]]) +
		(A[i][k] - A[j][k])*(B[p[j]][p[k]] - B[p[i]][p[k]]);
	return(d);
}

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
long compute_delta_part(type_vector & p, type_matrix & delta,
	int i, int j, int r, int s)
{
	return(delta[i][j] + (A[r][i] - A[r][j] + A[s][j] - A[s][i])*
		(B[p[s]][p[i]] - B[p[s]][p[j]] + B[p[r]][p[j]] - B[p[r]][p[i]]) +
		(A[i][r] - A[j][r] + A[j][s] - A[i][s])*
		(B[p[i]][p[s]] - B[p[j]][p[s]] + B[p[j]][p[r]] - B[p[i]][p[r]]));
}


void intializeFittness(int index)
{
  int current_cost = 0;
  for (int i = 0; i < variablesize; i = i + 1) 
	for (int j = 0; j < variablesize; j = j + 1)
    {
		current_cost = current_cost + A[i][j] * B[allswarm[index][i]][allswarm[index][j]];
		if (i < j) 
		{
			deltaArray[index][i][j] = compute_delta(allswarm[index], i, j);
		}
    }
  for(int j = 0 ; j< variablesize ; j++)
  {
	  allbaseswarm[index][j] = allswarm[index][j];
  }
  if(current_cost < allbestfittness[index])
	 {
		 for(int j = 0 ; j< variablesize ; j++)
		 {
			allbest[index][j] = allswarm[index][j];
		 }
		 allbestfittness[index] = current_cost;
		 if(current_cost < systembestfitness)
		 {
			 bestfittnessupdated = true;
			 systembestfitness = current_cost;
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresult[j] = allswarm[index][j];
			 }
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresultlist[systemIterator][j] = allswarm[index][j];
			 }
			 systemfunctioncounter[systemIterator] = functioncounter;
			 systemParticleIndex[systemIterator] = index;
			 systembestfitnesslist[systemIterator] = current_cost;
			 systemtimeinSec[systemIterator] = swatch.get_time_so_far("FullRun");
			 systemIterator ++; 
		 }
	 }
  allfittness[index] =  current_cost;
}

void intializeSwarms()
 {
	 systembestfitness = LONG_MAX;
	 for (int i = 0; i < particlesize; i++)
     {
			gererateRandom(i);
            for(int j = 0 ; j< variablesize ; j++)
			{
				allbest[i][j] = allswarm[i][j];
			}
			allbestfittness[i] = LONG_MAX;
            generateRandomVelocity(i);
            intializeFittness(i);
     }
 }

 int getParent(int i)
 {
            int result;
            if (allparticalparentindex[i] == -1)
            {
                result = 0;
            }
            else
            {
                result = allparticalparentindex[i];
            }
            return result;
 }

 void calculateVelocity(int particalIndex)
 {
	 
	 float rand1 = (float)ran01(&seed);
     float rand2 = (float)ran01(&seed);
     float rand3 = (float)ran01(&seed);
      
		
      //get the index of the parent particle
     int parentIndex = getParent(particalIndex);
     //multiple the current velocity with Inertia
	 for(int i=0;i< variablesize;i++)
	 {
		 for(int j=0;j<variablesize ; j++)
		 {
			 allvelocity[parentIndex][i][j] = allvelocity[parentIndex][i][j] * inertia; 
		 }
	 }
     // calculate the diffrence between current particle and its best result found yet
   double positionvalue = 0;
   double bestpositionvalue = 0;
   for(int i=0;i< variablesize;i++)
	 {
		 for(int j=0;j<variablesize ; j++)
		 {
			 positionvalue = 0;
			 bestpositionvalue = 0;
			 if(allswarm[particalIndex][i] == j)
				 positionvalue = 1;
			 if(allbest[particalIndex][i] == j)
				 bestpositionvalue = 1;
			 allvelocity[particalIndex][i][j] += (bestpositionvalue - positionvalue) * rand1 * cof1;
		 }
	 }
            // calculate the difference between the parent particle position anf the particle result
    
   double parentpositionvalue = 0;
   for(int i=0;i< variablesize;i++)
	 {
		 for(int j=0;j<variablesize ; j++)
		 {
			 positionvalue = 0;
			 parentpositionvalue = 0;
			 if(allswarm[particalIndex][i] == j)
				 positionvalue = 1;
			 if(allswarm[parentIndex][i] == j)
				 parentpositionvalue = 1;
			 allvelocity[particalIndex][i][j] += (parentpositionvalue - positionvalue) * rand1 * cof1;
		 }
	 }
   //clip velocity to Max and Min Values
      for(int i=0;i< variablesize;i++)
		 {
			 for(int j=0;j<variablesize ; j++)
			 {
				 if(allvelocity[particalIndex][i][j] <= minvelocity)
					 allvelocity[particalIndex][i][j] = minvelocity;
				 if(allvelocity[particalIndex][i][j] >= maxvelocity)
					 allvelocity[particalIndex][i][j] = maxvelocity;
			 }
		  }
 }

 void generate_random_solution(int* & p)
 {int i;
 for (i = 0; i < variablesize; i = i+1) p[i] = i;
 for (i = 0; i <  variablesize; i = i+1) transpose(p[i], p[unif(i,variablesize -1)]);
 }

 void recalculate(int particleIndex)
 {
	 int* solutionGenerate = new int [variablesize];
	 int* selectedresults = new int [variablesize];
	 for( int i=0 ; i< variablesize ; i++)
	 {
		 selectedresults[i] = 0;
	 }
	 generate_random_solution(solutionGenerate);
	 double* currentlist = new double [variablesize];
	 double result = 0;
	 double sumation = 0;
	 int count = 0 ;
	 int currentposition = 0;
	 double acumlatesum = 0;
	 double random = 0;
	 int k =0;
	 for( int i=0 ; i< variablesize ; i++)
	 {
		 sumation = 0;
		 count = 0;
		 
		 currentposition = solutionGenerate[i];
		 for(  k = 0 ; k < variablesize ; k++)
		 {
			 result = _I32_MIN;
			 if(selectedresults[k] == 0 )
			 {
				 result = allvelocity[particleIndex][currentposition][k];


				 result = 1.0 / ( 1.0 + exp(-result));
				
			 }
			  currentlist[k] = result;
		 }
		
		 //sumation of the velocity values of all particles
		 for ( k = 0; k < variablesize; k++)
         {
			 if(selectedresults[k] == 0)
			 {
				sumation += currentlist[k];
				count++;
			 }
         }

		 //sumation of the velocity values of all particles
		 for ( k = 0; k < variablesize; k++)
         {
			 if(selectedresults[k] == 0 )
			 {
				 currentlist[k] = currentlist[k]/ sumation;
			 }
         }

		 acumlatesum = 0;
	     random = ran01(&seed);
		for ( k = 0; k < variablesize; k++)
		{
			 if(selectedresults[k] == 0 )
			 {
				 acumlatesum += currentlist[k];
				 if(acumlatesum >= random)
				 {
					 allswarm[particleIndex][currentposition] = k;
					 selectedresults[k] = 1;
					 break;
				 }
			 }
		}	 
	 }
	
	 delete[] solutionGenerate;
	 delete[] currentlist;
	 delete[] selectedresults;
 }

 void swapParticals(int i, int p)
 {
      int* besttemp;
	  int* basetemp;
      int* swarmtemp;
      double** velocitytemp;
	  long** deltatemp;
      long fittnesstemp = 0;
      long bestfitnesstemp = 0;

      besttemp = allbest[i];
      swarmtemp = allswarm[i];
      velocitytemp = allvelocity[i];
      fittnesstemp = allfittness[i];
      bestfitnesstemp = allbestfittness[i];
	  deltatemp = deltaArray[i];
	  basetemp = allbaseswarm[i];


      allbestfittness[i] = allbestfittness[p];
      allbest[i] = allbest[p];
      allswarm[i] = allswarm[p];
      allvelocity[i] = allvelocity[p];
      allfittness[i] = allfittness[p];
	  deltaArray[i] = deltaArray[p];
	  allbaseswarm[i] = allbaseswarm[p];

      allbest[p] = besttemp;
      allswarm[p] = swarmtemp;
      allvelocity[p] = velocitytemp;
      allfittness[p] = fittnesstemp;
      allbestfittness[p] = bestfitnesstemp;
	  deltaArray[p] = deltatemp;
	  allbaseswarm[p]  = basetemp;
	  
 }

 void reArrangeBreadthFirst()
 {
	 int index = -1; 
	 for(int i=0;i < particlesize - degree ; i++)
	 {
		 long max = LONG_MAX;
		 for( int j=0; j < particlesize ; j++)
		 {
			  if (allparticalparentindex[j] == i)
			  {
				  if(max > allbestfittness[j])
				  {
					  index = j;
					  max = allbestfittness[j];
				  }
			  }
		 }
		  // if the best preforming children fittness better than its fittness swap particles
		if (allbestfittness[i] > allbestfittness[index])
		{
			 swapParticals(i, index);
		}
	 }
 }

void readResults(string filename)
{
	ifstream read;
	read.open(filename);
	int variable;
	read>>variable;
	read>>bestfitnessknown;
	read.close();
}

void buildSystem()
{
     intializeheirarchy();
     intializeSwarms();
}

void printout(int index, string word)
{
	/*cout<<word;
	cout<<endl;
	for(int i=0;i<variablesize;i++)
		cout<<setw(3)<<allswarm[index][i];
	cout<<endl;*/
}

long GACheckFitnessFor(int Ind)
{
  long total = 0;
  for (int i = 0; i < variablesize; i++)
  {
    for (int j = 0; j < variablesize; j++)
     {
		 total += (A[i][j]) * (B[Children[Ind][i]][Children[Ind][j]]); 
     }
  }
  return total;
}

void GASwapMutate(int Ind)
{
	int BitIndex1 = rand() % variablesize;
	int BitIndex2 = rand() % variablesize;

	int t = Children[Ind][BitIndex1];
	Children[Ind][BitIndex1] = Children[Ind][BitIndex2];
	Children[Ind][BitIndex2] = t;

}

void GAReverseMutate(int Ind)
{
	int BitIndex1 = rand() % variablesize;
	int BitIndex2 = rand() % variablesize;
	int interval = abs(BitIndex2 - BitIndex1) + 1;
	
	int* RevInterval = new int[interval];
	
	if(BitIndex1 > BitIndex2)
	{
		int t = BitIndex1;
		BitIndex1 = BitIndex2;
		BitIndex2 = t;
	}
	int count = BitIndex2;
	for(int i = 0; i <interval; i ++)
	{
		RevInterval[i] = Children[Ind][count];
		count --;
	}
	count = 0;
	for(int i = BitIndex1; i <= BitIndex2; i ++)
	{
		Children[Ind][i] = RevInterval[count];
		count++;
	}

	
	delete[] RevInterval;
}

void GACrossOver()
{
	int BitIndex = rand() % variablesize;

	/*int** NewCh = new int*[2];
	for(int i = 0; i < 2; i++)
		NewCh[i] = new int[variablesize];*/

	for( int i = 0; i < BitIndex; i++)
	{
		Children[0][i] = SelectedParents[0][i];
		Children[1][i] = SelectedParents[1][i];
	}
	int count1 = BitIndex;
	int count2 = BitIndex;
	for( int i = 0; i < variablesize; i++)
	{
		if(!ifContains(Children[1],SelectedParents[0][i], variablesize))
		{
			Children[1][count1] = SelectedParents[0][i];
			count1++;
		}
		if(!ifContains(Children[0],SelectedParents[1][i], variablesize))
		{
			Children[0][count2] = SelectedParents[1][i];
			count2++;
		}
	}

}

void GASortThisFit( int size)
{
	float T; int TT, n = size - 1;
	bool flag;
	do
	{ flag = false;
	for (int j = 0; j < n; j++)
		if (FitnessComp[j] > FitnessComp[j+1])
		{
			T = FitnessComp[j];
			FitnessComp[j] = FitnessComp[j+1];
			FitnessComp[j+1] = T;

			TT = Index[j];
			Index[j] = Index[j+1];
			Index[j+1] = TT;

			flag = true;
		}
	} while (flag);
}

void GASortThisSelecPar(int size)
{

	float T; int TT, n = size - 1;
	bool flag;
	do
	{ flag = false;
	for (int j = 0; j < n; j++)
		if (SelectedFitness[j] > SelectedFitness[j+1])
		{
			T = SelectedFitness[j];
			SelectedFitness[j] = SelectedFitness[j+1];
			SelectedFitness[j+1] = T;

			TT = SelectedIndexes[j];
			SelectedIndexes[j] = SelectedIndexes[j+1];
			SelectedIndexes[j+1] = TT;

			flag = true;
		}
	} while (flag);
}

void GASelectParents()
{
	//Tournament

	//int* ParentIndexes = new int[2];

	int Interval = unif(2, particlesize -1);

	int startingIndex = unif(0, particlesize-Interval);

	int endingIndex =  startingIndex + Interval;

	SelectedFitness = new long[Interval];
	SelectedIndexes = new int[Interval];
	ParentIndexes = new int[2];


	for(int i=0; i <Interval; i++)
	{
		SelectedIndexes[i] = i;
		SelectedFitness[i] = CurrentGenFitness[startingIndex + i];
	}
	GASortThisSelecPar(Interval);

	ParentIndexes[0] = startingIndex + SelectedIndexes[0];
	ParentIndexes[1] = startingIndex + SelectedIndexes[1];

	
	delete[] SelectedFitness;
	delete[] SelectedIndexes;
}

void GAStartPopulating()
{
	int newIndividuals = 0;
	float currentPm = ran01(&seed);
	float currentPc = ran01(&seed);
	
	bool newchildren;

	
	
	

	while(newIndividuals < particlesize)
	{
//		cout<<"Individual NO: "<<newIndividuals<<endl;
		newchildren = false;
		currentPm = ran01(&seed);
		currentPc = ran01(&seed);
		

		
		for(int i = 0; i < 2; i++)
			for(int j=0;j<variablesize;j++)
				Children[i][j] = -1;
		
		GASelectParents();
		for(int i = 0; i < variablesize ; i ++)
		{
			SelectedParents[0][i] = CurrentGeneration[ParentIndexes[0]][i];
			SelectedParents[1][i] = CurrentGeneration[ParentIndexes[1]][i];

		}
		if( currentPc <= CrossoverProb)
		{
			GACrossOver();
			newchildren = true;
		}
		else
			for(int i = 0; i < 2; i++)
				for(int j = 0; j < variablesize; j ++)
					Children[i][j] = SelectedParents[i][j];


		if( currentPm <= MutationProb)
		{
			int MutationChoice = unif(1,2);
			if(MutationChoice == 1)
				GASwapMutate(0);
			else
				GAReverseMutate(0);
			MutationChoice = unif(1,2);
			if(MutationChoice == 1)
				GASwapMutate(1);
			else
				GAReverseMutate(1);
			
			newchildren = true;
		}



		if(newchildren == false)
		{
			for(int j = 0; j < variablesize; j ++)
					NewPopulation[newIndividuals][j] = SelectedParents[0][j];
			NewPopFitness[newIndividuals] = CurrentGenFitness[ParentIndexes[0]];
			newIndividuals++;
			//NewPopulation[newIndividuals] = SelectedParents[1];
			//NewPopFitness[newIndividuals] = CurrentGenFitness[ParentIndexes[1]];
			//newIndividuals++;
		}
		else
		{
			FitnessComp[0] = GACheckFitnessFor(0);
			FitnessComp[1] = GACheckFitnessFor(1);
			FitnessComp[2] = CurrentGenFitness[ParentIndexes[0]];
			FitnessComp[3] = CurrentGenFitness[ParentIndexes[1]];

			
			Index[0] = -1;
			Index[1] = -2;
			Index[2] = ParentIndexes[0];
			Index[3] = ParentIndexes[1];
			
			//taking the best 1 out of the 4

			GASortThisFit(4);
			for(int i = 0; i < 1; i ++)
			{
				if(Index[i]>=0)			
				{
					for(int j = 0; j < variablesize; j ++)
						NewPopulation[newIndividuals][j] = CurrentGeneration[Index[i]][j];
					NewPopFitness[newIndividuals] = FitnessComp[i];
				}
				else if(Index[i] == -1)
				{
					for(int j = 0; j < variablesize; j ++)
						NewPopulation[newIndividuals][j] = Children[0][j];
					NewPopFitness[newIndividuals] = FitnessComp[i];
				}
				else if(Index[i] == -2)
				{
					for(int j = 0; j < variablesize; j ++)
						NewPopulation[newIndividuals][j] = Children[1][j];
					NewPopFitness[newIndividuals] = FitnessComp[i];
				}
				newIndividuals++;
			}
		}	
		
	}


}

void InitalizeGA()
{
	CurrentGeneration = new int*[particlesize];
	for(int i = 0; i < particlesize; i++)
		CurrentGeneration[i] = new int[variablesize];

	NewPopulation = new int*[particlesize];
	for(int i = 0; i < particlesize; i++)
		NewPopulation[i] = new int[variablesize];
	
	SelectedParents = new int*[2];

	for(int i = 0; i < 2; i++)
		SelectedParents[i] = new int[variablesize];
	Children = new int*[2];

	for(int i = 0; i < 2; i++)
		Children[i] = new int[variablesize];



	CurrentGenFitness = new long[particlesize];

	NewPopFitness = new long[particlesize];

	MutationProb = 0.3;
	CrossoverProb = 0.7;
	
}

void DisposeGA()
{

	for (int i = 0; i<particlesize; i++)
		delete[] CurrentGeneration[i];
	delete[] CurrentGeneration;
	
	for (int i = 0; i<particlesize; i++)
		delete[] NewPopulation[i];
	delete[] NewPopulation;

	for (int i = 0; i<2; i++)
		delete[] Children[i];
	delete[] Children;
	
	for (int i = 0; i<2; i++)
		delete[] SelectedParents[i];
	delete[] SelectedParents;


	delete[] CurrentGenFitness;

	delete[] NewPopFitness;

}

void GA()
{
	for(int i = 0; i < particlesize; i++)
	{
		CurrentGenFitness[i] = allbestfittness[i];
		for(int j= 0 ; j < variablesize; j ++)
			CurrentGeneration[i][j] = allbest[i][j];
	}
	GAStartPopulating();
	for(int i = 0; i < particlesize; i++)
	{
		allbestfittness[i] = NewPopFitness[i];
		if(NewPopFitness[i] < systembestfitness)
		{
			 bestfittnessupdated = true;
			 systembestfitness = NewPopFitness[i];
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresult[j] = NewPopulation[i][j];
			 }
			 for(int j = 0 ; j< variablesize ; j++)
			 {
				systembestresultlist[systemIterator][j] = NewPopulation[i][j];
			 }
			 systemfunctioncounter[systemIterator] = iteration;
			 systemParticleIndex[systemIterator] = -1;
			 systembestfitnesslist[systemIterator] = NewPopFitness[i];
			 systemtimeinSec[systemIterator] = swatch.get_time_so_far("FullRun");
			 systemIterator ++; 
		}
		for(int j= 0 ; j < variablesize; j ++)
			allbest[i][j] = NewPopulation[i][j];
	}
}

void GAParticle(int index)
{
	
	for(int i = 0; i < particlesize; i++)
	{
		CurrentGenFitness[i] = allbestfittness[i];
		for(int j= 0 ; j < variablesize; j ++)
			CurrentGeneration[i][j] = allbest[i][j];
	}
	double currentPm = ran01(&seed);
	double currentPc = ran01(&seed);
	bool newchildren = false;

	for(int i = 0; i < 2; i++)
			for(int j=0;j<variablesize;j++)
				Children[i][j] = -1;

	GASelectParents();
	for(int i = 0; i < variablesize ; i ++)
	{
		SelectedParents[0][i] = allswarm[index][i];
		SelectedParents[1][i] = CurrentGeneration[ParentIndexes[0]][i];
	}
	if( currentPc <= CrossoverProb)
	{
		GACrossOver();
		newchildren = true;
	}
	else
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < variablesize; j ++)
				Children[i][j] = SelectedParents[i][j];


	if( currentPm <= MutationProb)
	{
		int MutationChoice = unif(1,2);
		if(MutationChoice == 1)		
			GASwapMutate(0);
		else
			GAReverseMutate(0);
		MutationChoice = unif(1,2);
		if(MutationChoice == 1)
			GASwapMutate(1);
		else
			GAReverseMutate(1);
			
		newchildren = true;
	}
	if(newchildren == false)
	{
		for(int j = 0; j < variablesize; j ++)
			allswarm[index][j] = SelectedParents[0][j];
	}
	else
	{
			FitnessComp[0] = GACheckFitnessFor(0);
			FitnessComp[1] = GACheckFitnessFor(1);
			FitnessComp[2] = CurrentGenFitness[ParentIndexes[0]];
			FitnessComp[3] = CurrentGenFitness[ParentIndexes[1]];

			
			Index[0] = -1;
			Index[1] = -2;
			Index[2] = ParentIndexes[0];
			Index[3] = ParentIndexes[1];
			
			//taking the best 2 out of the 4

			GASortThisFit(4);
			for(int i = 0; i < 1; i ++)
			{
				if(Index[i]>=0)			
				{
					for(int j = 0; j < variablesize; j ++)
						allswarm[index][j] = CurrentGeneration[Index[i]][j];
				}
				else if(Index[i] == -1)
				{
					
					for(int j = 0; j < variablesize; j ++)
						allswarm[index][j] = Children[0][j];
				}
				else if(Index[i] == -2)
				{
					for(int j = 0; j < variablesize; j ++)
						allswarm[index][j] = Children[1][j];
				}
				
			}
	}
	

}

void Recalculateprobability()
{
	if(GAAlgorithm > HPSOAlgorithm)
	{
		if(probabilityGA < 0.75)
			probabilityGA += probabilityGA * 0.01;
	}
	if(GAAlgorithm < HPSOAlgorithm)
	{
		if(probabilityGA > 0.25)
			probabilityGA -= probabilityGA * 0.01;
	}
	PEAValues[iteration] = probabilityGA;
	GACount[iteration] = GAAlgorithm;
	PACount[iteration] = HPSOAlgorithm;
}

void Run()
 {
	 particalechosenGA = false;
	 GAAlgorithm = 0;
	 HPSOAlgorithm = 0;
	 if(HPSO)
	 {
		 for (int i = 0; i < particlesize; i++)
		 {
			 particalechosenGA = false;
			 if (hybridgeneticalgorithm)
			 {
				 double random = ran01(&seed);
				 if (random < probabilityGA)
					 particalechosenGA = true;
			 }
			 if (particalechosenGA)
			 {
				 GAParticle(i);
			 }
			 else
			 {
				 calculateVelocity(i);
				 recalculate(i);
			 }


			 calculateFittness(i, particalechosenGA);

			 if (systembestfitness <= bestfitnessknown)
				 break;
		 }
		
		 
		 reArrangeBreadthFirst();
	 }
	 if(geneticalgorithm)
		GA();
	 if(adaptivehybridgeneticalgorithm)
		 Recalculateprobability();
 }

void RunSystem()
{
	maxiterations = 100000 * variablesize;
	SizeData = maxiterations /4;
	iteration = 0;
	int failures = 0;
	int maxfailures = 100 * variablesize;
	intializeSwarm();
	probabilityGA = 0.7;
	GACount = new int[maxiterations];
	PACount = new int[maxiterations];
	PEAValues = new double[maxiterations];
	buildSystem();
	InitalizeGA();
	while(iteration < maxiterations)
	{
		bestfittnessupdated = false;
		//InitalizeGA();
		Run();

		if(!bestfittnessupdated)
			failures ++;
		else
			failures = 0;

		if(failures > maxfailures)
			break;

		if(systembestfitness <= bestfitnessknown)
			break;

		if(swatch.get_time_so_far("FullRun") > 10800)
			break;

		iteration++;

	}
	
}

string GetExePath() 
{
	char buffer[MAX_PATH];
    GetModuleFileNameA( NULL, buffer, MAX_PATH );
    string::size_type pos = string( buffer ).find_last_of( "\\/" );
    return string( buffer ).substr( 0, pos);
}

string algorithmName() 
{
	string algorithmName = "";
	if(geneticalgorithm || hybridgeneticalgorithm)
		algorithmName += "GA";
	if(HPSO)
		algorithmName +=  "HPSO";
	if(adaptivehybridgeneticalgorithm)
		algorithmName += " Adaptive";
	if(hybridgeneticalgorithm)
		algorithmName += " Hybrid";
	
	return algorithmName;
}

void main(int argc, char *argv[])
{
	string filename = argv[1];
	string resultsfilename =  argv[2];

    unsigned found = filename.find_last_of("/\\");

    string partfile_name = filename.substr(found+1);

    string problemname = partfile_name.substr(0,partfile_name.find("."));
	string folderLog;
	seed = atoi(argv[3]);
	
	
	geneticalgorithm = atoi(argv[4]);
	hybridgeneticalgorithm = atoi(argv[5]);
	adaptivehybridgeneticalgorithm =  atoi(argv[6]);
	HPSO =  atoi(argv[7]);
	
	int procedure =   atoi(argv[8]);
	
	seed += procedure;

	swatch.set_mode(CPU_TIME);

	readFile(filename);
	readResults(resultsfilename);
	
	swatch.start("FullRun");
	RunSystem();
	swatch.stop("FullRun");

	string folder = GetExePath() + "\\" + "Output";

	CreateDirectoryA(folder.c_str(), NULL);

	folder += "\\" + algorithmName();

	CreateDirectoryA(folder.c_str(), NULL);
	
	folderLog += "\\log";

	CreateDirectoryA(folderLog.c_str(), NULL);

	GUID guid;
	CoCreateGuid(&guid);

	stringstream outputfilestrem; 
	stringstream LogFileStream;

	outputfilestrem << folder;
	outputfilestrem << "\\";
	outputfilestrem << guid.Data1;
	outputfilestrem << guid.Data2;
	outputfilestrem << guid.Data3;
	for(int i = 0 ; i < 8 ; i++)
		outputfilestrem << (int)guid.Data4[i];
	outputfilestrem << "_" ;
	outputfilestrem << problemname;
	outputfilestrem <<".txt";

	/*LogFileStream << folderLog;
	LogFileStream << "\\";
	LogFileStream << guid.Data1;
	LogFileStream << guid.Data2;
	LogFileStream << guid.Data3;
	for(int i = 0 ; i < 8 ; i++)
		LogFileStream << (int)guid.Data4[i];
	LogFileStream << "_" ;
	LogFileStream << problemname;
	LogFileStream <<".log";*/
	
	ofstream myfile;
	myfile.open(outputfilestrem.str());

	//ofstream mylog;
	//mylog.open(LogFileStream.str());

	for( int i=0;i<systemIterator;i++)
	{
		myfile << setw(10)<< systembestfitnesslist[i]<< setw(10) << systemParticleIndex[i] << setw(10) << systemfunctioncounter[i] << setw(10)<< systemtimeinSec[i] << setw(5);
		for( int j= 0; j< variablesize;j++)
		{
			myfile<< setw(5) << systembestresultlist[i][j];
		}
		myfile<<endl;
	}

	//mylog << setw(10)<< "GA count"<< setw(10) <<"HPSO count" << setw(10) << "PEA" <<endl;
	//for( int i=0;i<iteration;i++)
	//{
	//	mylog << setw(10)<< GACount[i]<< setw(10) << PACount[i] << setw(10) << PEAValues[i] <<endl;
	//}
	//mylog<<endl;
	
	myfile.close();
//	mylog.close();
	Dispose();
	DisposeGA();
}