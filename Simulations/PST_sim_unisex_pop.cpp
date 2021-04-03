#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
#include <random>

//#include <libiomp/omp.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

int main(int argc, char * argv[]) {
    
	// Import parameters for simulation
    
	int N = atof(argv[1]);			// number of individuals
	double muB = atof(argv[2]);		// neutral locus mutation rate per N reproductive events
	
	double sampleRate = atof(argv[3]);	// Time interval between samples
	double minTime = atof(argv[4]);		// min time sampled (assumed relaxation time)
	double maxTime = atof(argv[5]);		// max time sampled
	double NK = atof(argv[6]);
	double sigmaB = atof(argv[7]);
	double c = atof(argv[8]);

	int K = N/(double)NK;

	std::ostringstream outfileString;
	ofstream outfile;

	outfileString << "/home/mapc-gwac20-user/Dropbox/CoalescentMatingTypes/FinalisedVersions/Simulations/UnisexPopData/sigB_"<< sigmaB <<"_c_"<< c <<"_muB_"<< muB <<"_k_"<< K <<"_N_"<< N <<"_mint_"<< minTime <<"_maxt_"<< maxTime <<"_dt_"<< sampleRate <<".txt";
	outfile.open( (outfileString.str()).c_str() );

	cout << "/home/mapc-gwac20-user/Dropbox/CoalescentMatingTypes/FinalisedVersions/Simulations/UnisexPopData/sigB_"<< sigmaB <<"_c_"<< c <<"_muB_"<< muB <<"_k_"<< K <<"_N_"<< N <<"_mint_"<< minTime <<"_maxt_"<< maxTime <<"_dt_"<< sampleRate <<".txt" << endl;

	double T[3];
	T[0] = N*(1-sigmaB)*(1-muB);
	T[1] = N*(1-sigmaB)*muB + N*sigmaB*muB;
	T[2] = N*sigmaB*(1-muB);

	double TSUM = T[0]+T[1]+T[2];
	//cout << TSUM << endl;

	int n[N][2];
	for(int i=0; i<N; i++)
	{
		n[i][1] = 0;
		n[i][0] = floor( i / (double)NK );
	}
	// define value for next mutated allele  
	int nextMut = 1;

	double t=0;
	double dT=0;
	double tau;

	int endProg = 0;

	srand ( time(NULL) );
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double r1;
	double r2;
	double sum;
	int z;

	while( endProg == 0 )
	{

		// Note; broad transition classes indepedent of n, therefore don't need to recalculate each event

		//generate random numbers
	        r1 = 0;
 	       	while(r1 == 0)
        	{
        	    r1 = dist(engine);
        	}
        	r2 = 0;
        	while(r2 == 0)
        	{
        	    r2 = dist(engine);
        	}

		dT = (1/TSUM)*(log(1/r1));

		sum = 0;
		z = -1;

		while(r2*TSUM >= sum)
		{
			z++;
		  	sum+=T[z];
		}

		// Note that the amount of times that we cycle through the data depends on dT (i.e. time until next mutation) but we dont want to output the updated state and time repeatedly

		while( (t+dT) > tau )
    		{
			if ( (t+dT) <= (maxTime) &&  t >= minTime   )
        		{
				outfile << tau <<" ";
	            		for(int i=0; i< N; i++)
	            		{
	                		outfile << n[i][1] <<" ";
	            		}
	            		outfile <<" "<< endl;
			}
	                else if( t > maxTime )
	                {
				endProg = 1;
	                }
        
	               tau += sampleRate;
		}

        	// END: OUTPUT DATA -------------------------------------------------------------------------------------
		
		t = t + dT;
		//cout << "event was z=" << z <<endl;
		
		// if asexual reproduction event
		if( z == 0 )
		{
			// pick individual to reproduce asexually (any individual in population) 
			int zi = N*dist(engine);
			// pick individual to die (any individual in population)
			int zj = N*dist(engine);
			// replace zj with zi.
			n[zj][1] = n[zi][1];
		}
		// if mutation event 
		else if( z==1 )
		{
			// pick individual to mutate (any individual in population) 
			int zi = N*dist(engine);
			n[zi][1] = nextMut;
			nextMut = nextMut + 1;

			// Periodically reduce the mutation counter
			if( nextMut > 5000)
			{
				int minMut = nextMut;
				for(int mi = 0; mi < N; mi++)
				{
					if(n[mi][1] < minMut)
					{
						minMut = n[mi][1];
					}
				}

				for(int mi = 0; mi < N; mi++)
				{

					n[mi][1] = n[mi][1] - minMut;
					
				}

				nextMut = nextMut - minMut;
			}

		}
		// if sexual reproduction event 
		else if( z==2 )
		{
			// pick focal individual (any individual in population) 
			int zi = N*dist(engine);
			// pick partner (any individual in population)
			int zj = N*dist(engine);
			int zk = N*dist(engine);

			// pick inheritance from either parent
			double zr = dist(engine);
			if( zr < 0.5)
			{
				n[zk][1] = n[zi][1];
			}
			else
			{
				n[zk][1] = n[zj][1];
			}
		}


		

	}
  	
	outfile.close();

	return 0;

}
