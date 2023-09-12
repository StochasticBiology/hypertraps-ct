// this is the workhorse code for the HyperTraPS-CT algorithm
// it takes command line arguments that dictate the data file(s), the structure of the inference run, and various parameters
// the output file is a set of samples from the posterior distribution inferred over hypercube parameters

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define RND drand48()

// maximum number of datapoints (just for memory allocation)
#define _MAXN 2000

// number of trajectories N_h, and frequencies of sampling for posteriors and for output
int BANK = 200;
int TMODULE = 100;

int _EVERYITERATION = 0;

// control output
int VERBOSE = 0;
int SPECTRUM_VERBOSE = 0;
int APM_VERBOSE = 0;

// impose limits on integer val to be between lo and hi
void limiti(int *val, int lo, int hi)
{
  if(*val < lo) *val = lo;
  if(*val > hi) *val = hi;
} 

// impose limits on double val to be between lo and hi
void limitf(double *val, int lo, int hi)
{
  if(*val < lo) *val = lo;
  if(*val > hi) *val = hi;
}

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

int mypow2(int r)
{
  int s = 1;
  int i;
  for(i = 1; i <= r; i++)
    s *= 2;
  return s;
}

void OutputTransitions(char *besttransstr, double *ntrans, int LEN)
{
  FILE *fp;
  int i, j, k;
  int statedec;
  int state[LEN];
  double rate, totrate;
  
  fp = fopen(besttransstr, "w");
  fprintf(fp, "From To Probability\n");
  
  for(i = 0; i < mypow2(LEN); i++)
    {
      statedec = i;
      for(j = LEN-1; j >= 0; j--)
	{
	  if(statedec >= mypow2(j))
	    {
	      state[LEN-1-j] = 1;
	      statedec -= mypow2(j);
	    }
	  else
	    state[LEN-1-j] = 0;
	}

      totrate = 0;
      for(j = 0; j < LEN; j++)
	{
	  /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
	  if(state[j] == 0)
	    {
	      rate = ntrans[j*LEN+j];
	      for(k = 0; k < LEN; k++)
		rate += state[k]*ntrans[k*LEN+j];
	      rate = exp(rate);
	      totrate += rate;
	    }
	}

      for(j = 0; j < LEN; j++)
	{
	  /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
	  if(state[j] == 0)
	    {
	      rate = ntrans[j*LEN+j];
	      for(k = 0; k < LEN; k++)
		rate += state[k]*ntrans[k*LEN+j];
	      rate = exp(rate);
	      fprintf(fp, "%i %i %e\n", i, i+mypow2(LEN-1-j), rate/totrate);
	    }
	}
    }
  fclose(fp);
}


void OutputStates(char *beststatesstr, double *ntrans, int LEN)
{
  FILE *fp;
  int i, j, k, a;
  int statedec;
  int src, dest;
  int state[LEN];
  double rate, totrate;
  int *active, *newactive;
  double *probs;
  int nactive, newnactive;
  int level;
  int found;
  
  fp = fopen(beststatesstr, "w");
  fprintf(fp, "State Probability\n");

  probs = (double*)malloc(sizeof(double)*mypow2(LEN));
  active = (int*)malloc(sizeof(int)*mypow2(LEN));
  newactive = (int*)malloc(sizeof(int)*mypow2(LEN));
  for(i = 0; i < mypow2(LEN); i++)
    probs[i] = 0;
  level = 0;
  
  probs[0] = 1;
  
  active[0] = 0;
  nactive = 1;
  
  while(nactive > 0)
    {
      newnactive = 0;
      /*      printf("%i active\n", nactive);
	      for(a = 0; a < nactive; a++)
	      printf("%i ", active[a]);
	      printf("\n\n"); */
	    
      for(a = 0; a < nactive; a++)
	{
	  src = active[a];
	  statedec = src;
	  for(j = LEN-1; j >= 0; j--)
	    {
	      if(statedec >= mypow2(j))
		{
		  state[LEN-1-j] = 1;
		  statedec -= mypow2(j);
		}
	      else
		state[LEN-1-j] = 0;
	    }

	  totrate = 0;
	  for(j = 0; j < LEN; j++)
	    {
	      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
	      if(state[j] == 0)
		{
		  rate = ntrans[j*LEN+j];
		  for(k = 0; k < LEN; k++)
		    rate += state[k]*ntrans[k*LEN+j];
		  rate = exp(rate);
		  totrate += rate;
		}
	    }

	  for(j = 0; j < LEN; j++)
	    {
	      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
	      if(state[j] == 0)
		{
		  dest = src+mypow2(LEN-1-j);
		  rate = ntrans[j*LEN+j];
		  for(k = 0; k < LEN; k++)
		    rate += state[k]*ntrans[k*LEN+j];
		  rate = exp(rate);
		  probs[dest] += probs[src] * rate/totrate;
		  //		  printf("%i: %i (from %i, %e): %e\n", level, dest, src, probs[src], probs[dest]);
		
		  found = 0;
		  for(k = 0; k < newnactive; k++)
		    {
		      if(newactive[k] == dest) { found = 1; break; }
		    }
		  if(found == 0)
		    newactive[newnactive++] = dest;
		}
	    }
	}
      for(a = 0; a < newnactive; a++)
	active[a] = newactive[a];
      nactive = newnactive;
      level++;
    }
  
  for(dest = 0; dest < mypow2(LEN); dest++)
    fprintf(fp, "%i %e\n", dest, probs[dest]);

  
  fclose(fp);
  free(active);
  free(newactive);
  free(probs);
}


// pick a new locus to change in state "state"; return it in "locus" and keep track of the on-course probability in "prob". "ntrans" is the transition matrix
void PickLocus(int *state, double *ntrans, int *targ, int *locus, double *prob, double *beta, int LEN)
{
  int i, j;
  double rate[LEN];
  double totrate, nobiastotrate;
  double cumsum[LEN];
  double r;

  nobiastotrate = 0;

  /* compute the rate of loss of gene i given the current genome -- without bias */
  for(i = 0; i < LEN; i++)
    {
      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
      if(state[i] == 0)
	{
	  rate[i] = ntrans[i*LEN+i];
	  for(j = 0; j < LEN; j++)
	    rate[i] += state[j]*ntrans[j*LEN+i];
	  rate[i] = exp(rate[i]);
	}
      else /* we've already lost this gene */
	rate[i] = 0;

      /* roulette wheel calculations as normal */
      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
      nobiastotrate += rate[i];
    }

  totrate = 0;

  /* compute the rate of loss of gene i given the current genome -- with bias */
  for(i = 0; i < LEN; i++)
    {
      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
      if(state[i] == 0 && targ[i] != 0)
	{
	  rate[i] = ntrans[i*LEN+i];
	  for(j = 0; j < LEN; j++)
	    rate[i] += state[j]*ntrans[j*LEN+i];
	  rate[i] = exp(rate[i]);
	}
      else /* we've already lost this gene OR WE DON'T WANT IT*/
	rate[i] = 0;

      /* roulette wheel calculations as normal */
      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
      totrate += rate[i];
    }

  /* normalised, additive rates -- is this sensible? */
  for(i = 0; i < LEN; i++)
    cumsum[i] /= totrate;

  r = RND;
  for(i = 0; i < LEN-1; i++)
    {
      if(cumsum[i] < r && cumsum[i+1] > r) { break; }
    }

  *locus = i;

  *prob = totrate/nobiastotrate;
  *beta = nobiastotrate;
}


// compute HyperTraPS probability of a transition from "startpos" to "targ" given transition matrix "P"
double LikelihoodMultiple(int *targ, double *P, int LEN, int *startpos, double tau1, double tau2)
{
  int *bank;
  int n0, n1;
  double *reject;
  int i, j, r;
  int locus;
  int attempt[LEN];
  double min;
  double mean;
  double *prodreject;
  double summand[LEN];
  int fail, score;
  int *hits;
  double totalsum;

  // new variables
  double u, prob_path, vi, betaci, nobiastotrate;
  double analyticI1, analyticI2;
  double sumI1, sumI2;
  int n;
  double tmprate;
  double *recbeta;
  // nobiastotrate is retain to match role in PickLocus but basically corresponds to -u
  int exitcount = 0;
  
  // allocate memory for BANK (N_h) trajectories
  bank = (int*)malloc(sizeof(int)*LEN*BANK);
  reject = (double*)malloc(sizeof(double)*BANK);
  hits = (int*)malloc(sizeof(int)*BANK);
  prodreject = (double*)malloc(sizeof(double)*BANK);
  recbeta = (double*)malloc(sizeof(double)*LEN*BANK);

  // initialise each trajectory at the start state; count 0s and 1s
  for(i = 0; i < LEN*BANK; i++)
    bank[i] = startpos[i%LEN]; 
  n0 = 0;
  for(i = 0; i < LEN; i++)
    n0 += startpos[i];

  n1 = 0;
  for(i = 0; i < LEN; i++)
    n1 += (targ[i] == 1 || targ[i] == 2);

  if(n0 > n1)
    {
      // the target comes before the source
      printf("Wrong ordering, or some other problem with input file. Remember that as of 2021, data file rows should be ordered ancestor then descendant!\n");
      exit(0);
    }

  mean = 1;
  totalsum = 0;

  // check we're not already there
  fail = 0;
  for(i = 0; i < LEN; i++)
    fail += (targ[i] != startpos[i]);
  if(fail == 0) totalsum = 1;

  for(i = 0; i < BANK; i++)
    prodreject[i] = 1;

  // loop through the number of evolutionary steps we need to make
  for(r = n0; r < n1; r++)
    {
      for(i = 0; i < BANK; i++)
	hits[i] = 0;

      // loop through each trajectory
      for(i = 0; i < BANK; i++)
	{
	  for(j = 0; j < LEN; j++)
	    attempt[j] = bank[LEN*i+j];
	  // pick the locus to change at this step, and record the probability that we stay on track to the target
	  PickLocus(&bank[LEN*i], P, targ, &locus, &reject[i], &recbeta[LEN*i + (r-n0)], LEN);
	  bank[LEN*i+locus] = 1;

	  fail = 0;
	  // count whether we're there or not
	  for(j = 0; j < LEN; j++)
	    {
	      if(bank[LEN*i+j] != targ[j] && targ[j] != 2) fail = 1;
	    }
	  hits[i] += (1-fail);

	}

      // keep track of total probability so far, and record if we're there
      summand[r] = 0;
      for(i = 0; i < BANK; i++)
	{
	  prodreject[i] *= reject[i];
	  summand[r] += prodreject[i]*hits[i];
	}
      summand[r] /= BANK;

      totalsum += summand[r];
    }

  // if we're not using continuous time, avoid this calculation and just return average path probability
  if(tau1 == 0 && tau2 == INFINITY)
    {
      if(n0 == n1) prob_path = 1;
      else
	{
	  prob_path = 0;
	  for(n = 0; n < BANK; n++)
	    {
	      prob_path += prodreject[n]*1./BANK;
	    }
	}
          
      free(bank);
      free(reject);
      free(hits);
      free(prodreject);
      free(recbeta);

      return prob_path;
    }

  // we're using continuous time
  // now, compute loss intensity from this state
  nobiastotrate = 0;
  for(i = 0; i < LEN; i++)
    {
      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
      if(targ[i] == 0)
	{
	  tmprate = P[i*LEN+i];
	  for(j = 0; j < LEN; j++)
	    tmprate += targ[j]*P[j*LEN+i];
	  tmprate = exp(tmprate);
	}
      else /* we've already lost this gene */
	tmprate = 0;

      nobiastotrate += tmprate;
    }
  u = -nobiastotrate;

  // now go through the recorded paths and compute vi, betaci
  analyticI1 = analyticI2 = 0;
  if(n0 != n1)
    {
      for(n = 0; n < BANK; n++)
	{
	  prob_path = prodreject[n]*1./BANK;
          sumI1 = sumI2 = 0;
	  for(i = 0; i < n1-n0; i++)
	    {
	      vi = 1;
	      for(j = 0; j < n1-n0; j++)
		{
		  vi *= recbeta[LEN*n+j];
		  if(j != i)
		    vi *= 1./(recbeta[LEN*n+j]-recbeta[LEN*n+i]);
		  if(SPECTRUM_VERBOSE)
		    printf("step %i recbeta %.4f vi %.4f\n", j, recbeta[LEN*n+j], vi);
		}
	      betaci = recbeta[LEN*n+i];

	      sumI1 += vi*(exp(tau1*u)-exp(-tau1*(betaci)))/(betaci + u);
	      sumI2 += vi/betaci*(exp(-betaci*tau1) - exp(-betaci*tau2));

	      if(SPECTRUM_VERBOSE)
		printf("stepx %i vi %.4f betaci %.4f u %.4f | sumI1 %.4f sumI2 %.4f\n (tau1 %f tau2 %f)\n", i, vi, betaci, u, sumI1, sumI2, tau1, tau2);
	    }

	  if(sumI1 < 0 || sumI2 < 0)
	    {
	      //	      printf("I got a negative value for I1 (%e) or I2 (%e), which shouldn't happen and suggests a lack of numerical convergence. This can happen with large numbers of features. I'm stopping to avoid unreliable posteriors; consider running without continuous time option.\n", sumI1, sumI2);
	      //exit(0);
	    }
	  
	  analyticI1 += (prob_path*sumI1);
	  analyticI2 += (prob_path*sumI2);
	  if(SPECTRUM_VERBOSE)
	    printf("prob_path %.4f sumI1 %.4f sumI2 %.4f | analyticI1 %.4f analyticI2 %.4f\n", prob_path, sumI1, sumI2, analyticI1, analyticI2);
	  exitcount++;
	  //	  if(exitcount == 3)
	  //exit(0);
	}
    }
  else
    {
      analyticI1 += exp(u*tau1); // just probability of dwelling at start
    }

  free(bank);
  free(reject);
  free(hits);
  free(prodreject);
  free(recbeta);

  /*  if(analyticI1+analyticI2 > 100)
      {
      printf("Exiting at line 283\n");
      exit(0);
      }*/
  
  return analyticI1+analyticI2;
}

// get total likelihood for a set of changes
double GetLikelihoodCoalescentChange(int *matrix, int len, int ntarg, double *ntrans, int *parents, double *tau1s, double *tau2s)
{
  double loglik, tloglik, tlik;
  int i, j;
  int multiple;
  int startpos[len];

  // initialise and start at one corner of the hypercube
  loglik = 0;
  for(i = 0; i < len; i++)
    startpos[i] = 0;

  // loop through each ancestor/descendant pair
  for(i = 0; i < ntarg/2; i++)
    {
      // output if desired
      if(VERBOSE)
	{
	  printf("Target %i: ", i);
	  for(j = 0; j < len; j++) printf("%i", matrix[2*i*len+len+j]);
	  printf(" parent is: " );
	  for(j = 0; j < len; j++) {  startpos[j] = matrix[2*i*len+j]; printf("%i", startpos[j]); }
	  printf("\n");
	}
      // initialise start position
      for(j = 0; j < len; j++)
	startpos[j] = (matrix[2*i*len+j]);
      // get log-likelihood contribution from this pair (transition) using HyperTraPS
      tlik = LikelihoodMultiple(&(matrix[2*i*len+len]), ntrans, len, startpos, tau1s[i], tau2s[i]);
      tloglik = log(tlik);
      if(tlik < 0)
	{
	  printf("Somehow I have a negative likelihood, suggesting a lack of numerical convergence. Terminating to avoid unreliable posteriors.\n");
	  //exit(0);
	}

      // output if required
      if(VERBOSE)
	printf("--- %i %f %f\n", i, exp(tloglik), tloglik);
      loglik += tloglik;
    }

  // return total log likelihood
  return loglik;
}

void GetGradients(int *matrix, int len, int ntarg, double *trans, int *parents, double *tau1s, double *tau2s, double *gradients)
{
  double lik, newlik;
  int i;
  
  lik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
  for(i = 0; i < len*len; i++)
    {
      trans[i] += 0.1;
      newlik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
      gradients[i] = (newlik-lik)/0.1;
      trans[i] -= 0.1;
    }
}

void helpandquit(int debug)
{
  printf("Options [defaults]:\n\n--obs file.txt\t\tobservations file [NA]\n--times file.txt\t(start) timings file for CT [NA]\n--endtimes file.txt\tend timings file for CT [NT]\n--seed N\t\trandom seed [0]\n--length N\t\tchain length (10^N) [3]\n--kernel N\t\tkernel index [5]\n--walkers N\t\tnumber of walker samplers for HyperTraPS [200]\n--losses \t\tconsider losses (not gains) [OFF]\n--apm \t\t\tauxiliary pseudo-marginal sampler [OFF]\n--sgd\t\t\tuse gradient descent [OFF]\n--sa\t\t\tuse simulated annealing [OFF]\n--label label\t\tset output file label [OBS FILE AND STATS OF RUN]\n--help\t\t\t[show this message]\n--debug\t\t\t[show this message and detailed debugging options]\n\n");
  if(debug)
    printf("debugging options:\n--verbose\t\tgeneral verbose output [OFF]\n--spectrumverbose\tverbose output for CT calculations [OFF]\n--apmverbose\t\tverbose output for APM approach [OFF]\n--outputperiod N\tperiod of stdout output [100]\n(note: an undocumented option exists to pass CSV data as the observations file: file should have a header, and a column of (ignored) sample IDs, before subsequent columns with all \"before\" features followed by all \"after\" features on the same row.  \n\n");
  exit(0);
}

// main function processes command-line arguments and run the inference loop
int main(int argc, char *argv[])
{
  int parents[_MAXN];
  FILE *fp;
  int *matrix;
  int len, ntarg;
  double *trans, *ntrans, *gradients;
  int t;
  int i, j;
  char ch;
  double lik, nlik;
  int *rec, *tmprec;
  int maxt, allruns;
  int seed;
  char str[200];
  char fstr[200];
  char shotstr[200], bestshotstr[200], besttransstr[200], beststatesstr[200];
  double DELTA, MU;
  int NVAL;
  int expt;
  double acc, rej, lacc, lrej;
  int chain1, chain2;
  double prob;
  double *tmpmat;
  double r;
  char fstr1[100], fstr2[100];
  time_t timer;
  char buffer[25];
  struct tm* tm_info;
  double taus[_MAXN], tau1s[_MAXN], tau2s[_MAXN];
  int ntau;
  int nancount = 0;
  int spectrumtype;
  double bestlik = 0;
  int lengthindex, kernelindex;
  int SAMPLE;
  int losses;
  int apm_seed, old_apm_seed, apm_step;
  int apm_type;
  int csv;
  char likstr[100];
  double testval;
  char header[10000];
  char obsfile[1000], timefile[1000], endtimefile[1000];
  int searchmethod;
  int filelabel;
  char labelstr[1000];
  int crosssectional;
  int tmprow[1000];
  time_t start_t, end_t;
  double diff_t;
  struct timeval t_stop, t_start;
  
  printf("\nHyperTraPS(-CT)\nSep 2023\n\nUnpublished code -- please do not circulate!\nPublished version available at:\n    https://github.com/StochasticBiology/HyperTraPS\nwith stripped-down version at:\n    https://github.com/StochasticBiology/hypertraps-simple\n\n");

  // default values
  spectrumtype = 0;
  lengthindex = 3;
  kernelindex = 5;
  losses = 0;
  apm_type = 0;
  filelabel = 0;
  crosssectional = 0;
  searchmethod = 0;
  crosssectional = 0;
  strcpy(obsfile, "");
  strcpy(timefile, "");
  strcpy(endtimefile, "");

  // deal with command-line arguments
  if(argc < 2) helpandquit(0);
  for(i = 1; i < argc; i+=2)
    {
      if(strcmp(argv[i], "--obs\0") == 0) strcpy(obsfile, argv[i+1]);
      else if(strcmp(argv[i], "--label\0") == 0) { filelabel = 1; strcpy(labelstr, argv[i+1]); }
      else if(strcmp(argv[i], "--times\0") == 0) { spectrumtype = 1; strcpy(timefile, argv[i+1]); }
      else if(strcmp(argv[i], "--endtimes\0") == 0) strcpy(endtimefile, argv[i+1]);
      else if(strcmp(argv[i], "--seed\0") == 0) seed = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--length\0") == 0) lengthindex = atof(argv[i+1]);
      else if(strcmp(argv[i], "--kernel\0") == 0) kernelindex = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--losses\0") == 0) { losses = 1; i--;} 
      else if(strcmp(argv[i], "--apm\0") == 0) {apm_type = 1; i--;}
      else if(strcmp(argv[i], "--help\0") == 0) helpandquit(0);
      else if(strcmp(argv[i], "--debug\0") == 0) helpandquit(1);
      else if(strcmp(argv[i], "--verbose\0") == 0) { VERBOSE = 1; i--; }
      else if(strcmp(argv[i], "--crosssectional\0") == 0) { crosssectional = 1; i--; }
      else if(strcmp(argv[i], "--spectrumverbose\0") == 0) { SPECTRUM_VERBOSE = 1; i--; }
      else if(strcmp(argv[i], "--apmverbose\0") == 0) { APM_VERBOSE = 1; i--; }
      else if(strcmp(argv[i], "--outputperiod\0") == 0) TMODULE = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--walkers\0") == 0) BANK = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--sgd\0") == 0) { searchmethod = 1; i--; }
      else if(strcmp(argv[i], "--sa\0") == 0) { searchmethod = 2; i--; }
      
      else printf("Didn't understand argument %s\n", argv[i]);
    }
  limiti(&lengthindex, 0, 7);
  limiti(&kernelindex, 0, 7);

 
  if(spectrumtype == 1) {
    printf("Running HyperTraPS-CT with:\n[observations-file]: %s\n[start-timings-file]: %s\n[end-timings-file]: %s\n[random number seed]: %i\n[length index]: %i\n[kernel index]: %i\n[walkers]: %i\n[losses (1) or gains (0)]: %i\n[APM]: %i\n", obsfile, timefile, endtimefile, seed, lengthindex, kernelindex, BANK, losses, apm_type);
  } else {
    printf("Running HyperTraPS with:\n[observations-file]: %s\n[random number seed]: %i\n[length index]: %i\n[kernel index]: %i\n[walkers]: %i\n[losses (1) or gains (0)]: %i\n[APM]: %i\n", obsfile, seed, lengthindex, kernelindex, BANK, losses, apm_type);
  }
  switch(searchmethod) {
  case 0: printf("Using MH MCMC\n"); break;
  case 1: printf("Using SGD\n"); break;
  case 2: printf("Using SA\n"); break;
  } 
  
  // initialise and allocate
  maxt = pow(10, lengthindex);
  if(maxt <= 10000) SAMPLE = 100; else SAMPLE = 1000;

  if(_EVERYITERATION)
    SAMPLE = 1;

  srand48(seed);
  matrix = (int*)malloc(sizeof(int)*100000);

  // choose parameterisation based on command line
  expt = kernelindex;
  switch(expt)
    {
    case 0: DELTA = 0; break;
    case 1: DELTA = 0.005; MU = 0.1; break;
    case 2: DELTA = 0.05; MU = 1.; break;
    case 3: DELTA = 0.05; MU = 1.; break;
    case 4: DELTA = 0.1; MU = 1.; break;
    case 5: DELTA = 0.25; MU = 1.; break;
    case 6: DELTA = 0.5; MU = 1.; break;
    default: DELTA = 0.75; MU = 1.; break;
    }
  
  // read data on changes from input file
  // if we're thinking about losses, we're regarding gene losses as feature acquisitions; and thus inverting the data
  fp = fopen(obsfile, "r");
  if(fp == NULL)
    {
      printf("Couldn't find observations file %s\n", obsfile);
      return 0;
    }
  i = 0; len = 0; csv = 0;
  do{
    ch = fgetc(fp);
    if((ch != '0' && ch != '1' && ch != '2' && ch != ' ' && ch != '\t' && ch != '\n') && i == 0)
      {
	printf("Found non-digit character before any entries: interpreting as CSV file format\n");
	csv = 1;
	rewind(fp);
	do{ch = fgetc(fp); if(ch != '\n') header[i++] = ch; }while(ch != '\n');
	i = 0;
	ch = '\n';
      }
    switch(ch)
      {
      case '0': matrix[i++] = (losses == 1 ? 1 : 0); break;
      case '1': matrix[i++] = (losses == 1 ? 0 : 1); break;
      case '2': matrix[i++] = 2; break;
      case '\n':
	if(len == 0) len = i;
	if(csv) {
	  do{ch=fgetc(fp);}while(!feof(fp) && ch != ',');
	  do{ch=fgetc(fp);}while(!feof(fp) && ch != ',');
	}
	if(crosssectional) {
	  for(j = 0; j < len; j++)
	    {
	      tmprow[j] = matrix[i-len+j];
	      matrix[i-len+j] = 0;
	      matrix[i+j] = tmprow[j];
	    }
	  i += len;
	}
	break;
      }
  }while(!feof(fp));
  if(csv) len /= 2;
  ntarg = i/len;
  NVAL = len*len;
  fclose(fp);

  printf("Observed transitions:\n");
  for(i = 0; i < ntarg/2; i++)
    {
      for(j = 0; j < len; j++) printf("%i", matrix[2*len*i+j]);
      printf(" -> ");
      for(j = 0; j < len; j++) printf("%i", matrix[2*len*i+len+j]);
      printf("\n");
    }
  if(losses == 1) printf("(where 1 is absence)\n\n");
  if(losses == 0) printf("(where 1 is presence)\n\n");

  // grab timings from data file provided
  if(spectrumtype != 0)
    {
      ntau = 0;
      fp = fopen(timefile, "r");
      if(fp == NULL)
	{
	  printf("Couldn't find start timings file %s\n", timefile);
	  return 0;
	}
      while(!feof(fp))
	{
	  fscanf(fp, "%lf", &(tau1s[ntau]));
	  if(!feof(fp)) { ntau++; }
	}
      fclose(fp);
      if(ntau != ntarg/2) 
	{
	  printf("I found %i timings and %i observation pairs -- these numbers should be equal.\n", ntau, ntarg);
	  return 0;
	}

      fp = fopen(endtimefile, "r");
      if(fp == NULL)
	{
	  printf("Couldn't find end timings file -- I'm assuming that start times *are* end times (i.e. each observation has a precisely specified single time)\n");
	  for(i = 0; i < ntau; i++) tau2s[i] = tau1s[i];
	}
      else
	{
	  ntau = 0;
	  while(!feof(fp))
	    {
	      fscanf(fp, "%lf", &(tau2s[ntau]));
	      if(tau2s[ntau] < tau1s[ntau])
		{
		  printf("End time %f was less than start time %f!\n", tau2s[ntau], tau1s[ntau]);
		  exit(0);
		}
	      if(!feof(fp)) { ntau++; }
	    }
	  fclose(fp);
	}
      if(ntau != ntarg/2) 
	{
	  printf("I found %i timings and %i observation pairs -- these numbers should be equal.\n", ntau, ntarg);
	  return 0;
	}

    }
  else
    {
      ntau = ntarg/2;
      for(i = 0; i < ntau; i++)
	{
	  tau1s[i] = 0;
	  tau2s[i] = INFINITY;
	}
    }
  
  if(spectrumtype == 0)
    {
      printf("Number of features is %i, I found %i observation pairs\n", len, ntarg/2);
    }
  else
    {
      printf("Number of features is %i, I found %i observation pairs and %i timing pairs\n", len, ntarg/2, ntau);
      if(len > 30)
	{
	  printf("*** CAUTION: continuous time calculations sometimes fail to converge for large (>30) feature sets. This can lead to NaNs appearing, which will stop the simulation. Consider running without continuous time option.\n");
	}
    }
  printf("\n");
  
  // allocate memory and initialise output file
  trans = (double*)malloc(sizeof(double)*NVAL); 
  ntrans = (double*)malloc(sizeof(double)*NVAL);
  gradients = (double*)malloc(sizeof(double)*NVAL);
  tmpmat = (double*)malloc(sizeof(double)*NVAL);

  if(filelabel == 0)
    {
      sprintf(labelstr, "%s-%i-%i-%i-%i-%i-%i-%i", obsfile, spectrumtype, searchmethod, seed, lengthindex, kernelindex, BANK, apm_type);
    }
  // prepare output files
  sprintf(shotstr, "%s-posterior.txt", labelstr);
  fp = fopen(shotstr, "w"); fclose(fp);
  sprintf(bestshotstr, "%s-best.txt", labelstr);
  fp = fopen(bestshotstr, "w"); fclose(fp);
  sprintf(likstr, "%s-lik.txt", labelstr);
  fp = fopen(likstr, "w"); fprintf(fp, "Step,LogLikelihood1,LogLikelihood2\n"); fclose(fp);

  sprintf(besttransstr, "%s-trans.txt", labelstr);
  sprintf(beststatesstr, "%s-states.txt", labelstr);
  
  // initialise with an agnostic transition matrix
  for(i = 0; i < len*len; i++)
    trans[i] = 0;
  for(i = 0; i < len; i++)
    trans[i*len+i] = 1;

  // compute initial likelihood given this matrix
  time(&start_t);
  gettimeofday(&t_start, NULL);
  lik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
  time(&end_t);
  gettimeofday(&t_stop, NULL);
  diff_t = (t_stop.tv_sec - t_start.tv_sec) + (t_stop.tv_usec-t_start.tv_usec)/1.e6;
  //  diff_t = difftime(end_t, start_t);
  printf("One likelihood estimation took %e seconds.\n", diff_t);
        // MCMC or simulated annealing
      if(searchmethod == 0 || searchmethod == 2)
	{
	  printf("This code (%i steps) will probably take around %.3f seconds (%.3f hours) to complete.\n\n", maxt, diff_t*maxt, diff_t*maxt/3600.);
	}
      
  // initialise counters for acceptance ratio
  acc = rej = 0;
  lacc = lrej = 0;

  if(apm_type == 1)
    apm_seed = seed;
  
  // run the chain
  for(t = 0; t < maxt; t++)
    {
      if(t % SAMPLE == 0)
	{
	  // periodically output progress to a tracker file
	  time(&timer);
	  tm_info = localtime(&timer);

	  strftime(buffer, 25, "%Y:%m:%d %H:%M:%S", tm_info);

	  fp = fopen("alltrackernew.txt", "a");
	  fprintf(fp, "%s %i %i %i %s\n", shotstr, t, maxt/5, maxt, buffer);
	  fclose(fp);
	}
      // if we've got a new best likelihood, store it
      if(lik > bestlik || t == 0)
	{
	  bestlik = lik;
	  fp = fopen(bestshotstr, "w");
	  for(i = 0; i < len*len; i++)
	    fprintf(fp, "%f ", trans[i]);
	  fprintf(fp, "\n");
	  fclose(fp);

	  OutputTransitions(besttransstr, trans, len);
	  OutputStates(beststatesstr, trans, len);
	}

      // output some info periodically
      if(t % SAMPLE == 0)
	printf("%i - ", t);

      if(t > maxt/5 && t % SAMPLE == 0)
	{
	  // if we're burnt in, periodically sample the current parameterisation to an output file
	  // most appropriate for Bayesian MCMC but useful for all
	  fp = fopen(shotstr, "a");
	  for(i = 0; i < len*len; i++)
	    fprintf(fp, "%f ", trans[i]);
	  fprintf(fp, "\n");
	  fclose(fp);
	  fp = fopen(likstr, "a");
	  nlik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
	  fprintf(fp, "%i,%f,", t, nlik);
	  nlik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
	  fprintf(fp, "%f\n", nlik);
	  fclose(fp);
	}

      // MCMC or simulated annealing
      if(searchmethod == 0 || searchmethod == 2)
	{

	  if(apm_type == 0 || t%2 == 0)
	    {
	      // apply a perturbation to the existing parameterisation
	      // non-uniform priors can be employed here if desired 
		for(i = 0; i < NVAL; i++)
		  {
		    ntrans[i] = trans[i];
		    r = RND;
		    if(r < MU)
		      {
			ntrans[i] += gsl_ran_gaussian(DELTA);
		      }
		    if(ntrans[i] < -10) ntrans[i] = -10;
		    if(ntrans[i] > 10) ntrans[i] = 10;
		  }
	      if(APM_VERBOSE)
		{
		  printf("step 0 (change theta): apm_seed %i, ntrans[0] %f\n", apm_seed, ntrans[0]);
		}
	    }
	  else
	    {
	      // change the random number seed and keep the parameterisation the same
	      old_apm_seed = apm_seed;
	      apm_seed = seed+t;
	      for(i = 0; i < NVAL; i++)
		ntrans[i] = trans[i];
	      if(APM_VERBOSE)
		{
		  printf("step 1 (change u): apm_seed %i, ntrans[0] %f\n", apm_seed, ntrans[0]);
		}
	    }
      
	  // compute likelihood for the new parameterisation
	  if(apm_type == 1)
	    {
	      srand48(apm_seed);
	      if(APM_VERBOSE)
		{
		  printf("r seeded with %i, first call is %f\n", apm_seed, RND);
		}
	    }
	  nlik = GetLikelihoodCoalescentChange(matrix, len, ntarg, ntrans, parents, tau1s, tau2s);

	  if(APM_VERBOSE)
	    {
	      printf("likelihood %f\n", nlik);
	    }
      
	  // keep track of NaNs in calculations
	  if(isnan(nlik))
	    {
	      nancount++;
	    }

	  testval = RND;
	  if(searchmethod == 2)
	    {
	      testval = 0.1*sqrt(sqrt(t));
	    }

	  // compare likelihood to previous
	  if(nlik >= lik || -(lik-nlik) > log(testval))
	    {
	      // accept this new parameterisation
	      lik = nlik;
	  
	      if(apm_type == 0 || t%2 == 0)
		{
		  acc++; lacc++;
		  for(i = 0; i < NVAL; i++)
		    trans[i] = ntrans[i];
		}
	      if(APM_VERBOSE)
		{
		  printf("accepted: apm_seed %i trans[0] %f\n\n", apm_seed, trans[0]);
		}
	    }
	  else 
	    {
	      // reject the change
	      if(apm_type == 1 && t%2 == 1)
		{
		  apm_seed = old_apm_seed;
		}
	      else
		{
		  rej++; lrej++;
		}
	      if(APM_VERBOSE)
		{
		  printf("rejected: apm_seed %i trans[0] %f\n\n", apm_seed, trans[0]);
		}
	    }
	}
      // gradient descent
      if(searchmethod == 1)
	{
	  GetGradients(matrix, len, ntarg, trans, parents, tau1s, tau2s, gradients);
	  for(i = 0; i < len*len; i++)
	    {
	      trans[i] = trans[i]+gradients[i]*0.1;
	      if(trans[i] < -10) trans[i] = -10;
	      if(trans[i] > 10) trans[i] = 10;
	    }
	  
	  lik = GetLikelihoodCoalescentChange(matrix, len, ntarg, trans, parents, tau1s, tau2s);
	}
      //      if(t % SAMPLE == 0) printf("NaN count %i of %i\n", nancount, t);

      // output information periodically
      if(t % TMODULE == 0)
	{
	  printf("Iteration %i likelihood %f total-acceptance %f recent-acceptance %f trial-likelihood %f\n", t, lik, acc/(acc+rej), lacc/(lacc+lrej), nlik);
	  lacc = lrej = 0;
	}
    }

  return 0;
}
 

