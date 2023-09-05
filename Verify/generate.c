// generate synthetic datasets for testing HyperTraPS-CT code
// many time series simulated and sampled
// three different transition matrix choices: independent identical transition rates, and different modes of inter-trait influence

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

// number of traits
#define L 8

// get rates of transitions from the current state, given a particular choice of transition matrix
void ratesets(double *lambda, int expt, int *state)
{
  int i;
  
  switch(expt)
    {
      // independent, identical rates of trait acquisition
    case 0:
      for(i = 0; i < L; i++)
	lambda[i] = 0;
      break;

      // single path with some influences?
    case 1:
      lambda[0] = 1;
      for(i = 1; i <= 3; i++)
	lambda[i] = (state[i-1] == 1 ? 1 : -10);
      for(i = 4; i <= 5; i++)
	lambda[i] = (state[i-4] == 1 ? 1 : 0);
      for(i = 6; i < L; i++)
	lambda[i] = (state[i-6] == 1 ? (state[i-5] == 1 ? 1 : 0) : -10);
      break;

      // double path with some influences?
    case 2:
      lambda[0] = lambda[2] = 1;
      lambda[1] = (state[0] == 1 ? 1 : -10);
      lambda[3] = (state[2] == 1 ? 1 : -10);
      for(i = 4; i <= 5; i++)
	lambda[i] = (state[i-4] == 1 ? 1 : 0);
      for(i = 6; i < L; i++)
	lambda[i] = (state[i-6] == 1 ? (state[i-5] == 1 ? 1 : 0) : -10);
      break;
    }

  // make sure we can't acquire things we already have
  for(i = 0; i < L; i++)
    lambda[i] = exp(lambda[i])*(state[i] == 1 ? 0 : 1);
}
  
int main(int argc, char *argv[])
{
  int expt;
  int samples;
  int state[L], prevstate[L];
  double sampletimes[5];
  double t, maxt, dt;
  double cumsum[L], ratesum;
  double lambda[L];
  int i, j;
  int ball;
  double r;
  int nsamples;
  int ntimes;
  int change;
  double tmp;
  FILE *fp1, *fp2;
  char fstr[100];
  int featured[L];

  // sampling parameters
  nsamples = 1000;
  maxt = 1;
  ntimes = 1;

  // count the number of times each trait is featured
  for(i = 0; i < L; i++)
    featured[i] = 0;

  // allow random number seeding for reproducibility
  if(argc > 1)
    srand48(atoi(argv[1]));

  // loop through transition matrix choices
  for(expt = 0; expt <= 2; expt++)
    {
      // prep output files
      sprintf(fstr, "synth-%i-data.txt", expt);
      fp1 = fopen(fstr, "w");
      sprintf(fstr, "synth-%i-times.txt", expt);
      fp2 = fopen(fstr, "w");

      // run nsamples independent trajectories
      for(samples = 0; samples < nsamples; samples++)
	{
	  // initialise state
	  for(i = 0; i < L; i++)
	    state[i] = prevstate[i] = 0;

	  // choose sampling times for this trajectory, and bubble sort to order them
	  for(i = 0; i < ntimes; i++)
	    sampletimes[i] = RND*maxt;
	  do{
	    change = 0;
	    for(i = 0; i < ntimes-1; i++)
	      {
		if(sampletimes[i] > sampletimes[i+1])
		  {
		    tmp = sampletimes[i+1];
		    sampletimes[i+1] = sampletimes[i];
		    sampletimes[i] = tmp;
		    change = 1;
		  }
	      }
	  }while(change == 1);

	  // run through time in this trajectory
	  for(t = 0; t < maxt; )
	    {
	      // compute transition rates for this state
	      ratesets(lambda, expt, state);

	      // pick a transition
	      // build list of transitions and propensities
	      ratesum = 0; cumsum[0] = lambda[0];
	      for(i = 0; i < L; i++)
		{
		  if(i > 0)
		    cumsum[i] = cumsum[i-1]+lambda[i];
		  ratesum += lambda[i];
		}
	      // roll roulette ball and pick transition
	      r = RND*ratesum;
	      for(ball = 0; cumsum[ball] < r; ball++);

	      // time increment
	      dt = -log(RND) / ratesum;

	      // output if we've passed a sampling point
	      for(i = 0; i < ntimes; i++)
		{
		  // if this step will take us through a sampling point
		  if(t < sampletimes[i] && t+dt > sampletimes[i])
		    {
		      // screen output
		      printf("%i %i %f ", expt, samples, sampletimes[i]);
		      for(j = 0; j < L; j++)
			{
			  printf("%i", state[j]);
			}
		      printf("\n");

		      // output previous state to file, counting trait occurrences
		      for(j = 0; j < L; j++)
			{
			  fprintf(fp1, "%i ", prevstate[j]);
			  featured[j] += prevstate[j];
			}
		      fprintf(fp1, "\n");
		      // output current state to file, counting trait occurrences
		      for(j = 0; j < L; j++)
			{
			  fprintf(fp1, "%i ", state[j]);
      			  featured[j] += state[j];
			}
		      fprintf(fp1, "\n");

		      // output time to file
		      fprintf(fp2, "%f\n", (i == 0 ? sampletimes[i] : sampletimes[i]-sampletimes[i-1]));

		      // update previous state
		      for(j = 0; j < L; j++)
			prevstate[j] = state[j];
		    }
		}
	      // increment time and update state
	      t += dt;
	      state[ball] = 1;
	    }
	}
      fclose(fp1);
      fclose(fp2);
    }

  // report how many times each trait has featured as acquired (to ensure we don't get a randomly odd sampling set)
  int min = 9999999;
  int max = 0;
  for(i = 0; i < L; i++)
    {
      if(featured[i] < min) min = featured[i];
      if(featured[i] > max) max = featured[i];
      printf("%i ", featured[i]);
    }
  printf("\n%f\n", (double)(max-min)/min);
     
  return 0;
}
		  
    
  
