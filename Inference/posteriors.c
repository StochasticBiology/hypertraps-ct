// posteriors.c
// analyses samples from the posteriors produced by hypertraps-ct code and summarises ordering and time distributions
// note time distributions are only meaningful when the inference was performed using the continuous time options!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RND drand48()

// number of trajectories to simulate for each parameterisation
int NTRAJ = 100;

// number of times to call GetRoutes. each call runs NSAMP trajectories and records one sampled route, so the balance of the two controls the number of explicit recorded routes vs the number of sampled trajectories.
int NSAMP = 10;

// maximum continuous-time value above which results are truncated
#define MAXCT 1000

// just used in assigning ordinal labels to different features
#define FLEN 15

int mypow2(int r)
{
  int s = 1;
  int i;
  for(i = 1; i <= r; i++)
    s *= 2;
  return s;
}

int BinToDec(int *state, int LEN)
{
  int v = 1;
  int i;
  int val = 0;
  
  for(i = LEN-1; i >= 0; i--)
    {
      val += state[i]*v;
      v *= 2;
    }
  return val;
}

int nparams(int model, int LEN)
{
  switch(model)
    {
    case 0: return 0;
    case 1: return LEN;
    case 2: return LEN*LEN;
    case 3: return LEN*LEN*LEN;
    case 4: return LEN*LEN*LEN*LEN;
    case -1: return mypow2(LEN);
    default: return 0;
    }
}

double RetrieveEdge(int *state, int locus, double *ntrans, int LEN, int model)
{
  double rate;
  int i, j, k;
  
  if(model == 0)
    return 0;
  if(model == 1) // pi[locus] = rate of locus
    return ntrans[locus];
  if(model == 2) // pi[i*LEN + locus] = influence of i on locus
    {
      	      rate = ntrans[locus*LEN+locus];
	      for(i = 0; i < LEN; i++)
		rate += state[i]*ntrans[i*LEN+locus];
    }
  if(model == 3) // pi[j*LEN*LEN + i*LEN + locus] = influence of ij on locus, j>=i (j==i is influence of i on locus)
    {
      rate = ntrans[locus*LEN*LEN+locus*LEN+locus];
      for(i = 0; i < LEN; i++)
	{
	for(j = i; j < LEN; j++)
	  {
	    rate += state[i]*state[j]*ntrans[j*LEN*LEN+i*LEN+locus];
	  }
	}
    }
  if(model == 4) // pi[k*LEN*LEN*LEN + j*LEN*LEN + i*LEN + i] = influence of ijk on i, k>=j>=i (j==i is influence of ik, k==j==i is influence of i); what does k==j,i mean
    {
      rate = ntrans[locus*LEN*LEN*LEN+locus*LEN*LEN+locus*LEN+locus];
      for(i = 0; i < LEN; i++)
	{
	for(j = i; j < LEN; j++)
	  {
	    for(k = j; k < LEN; k++)
	      {
	    rate += state[i]*state[j]*state[k]*ntrans[k*LEN*LEN*LEN+j*LEN*LEN+i*LEN+locus];
	  }
	}
    }
    }
  if(model == -1)
    {
      rate = ntrans[BinToDec(state, LEN)*LEN+locus];
    }
  return exp(rate);
}

// simulate trajectories on a given hypercube parameterisation, and store a bunch of summary data about those trajectories
// mean[i] stores the mean acquisition ordering for feature i
// ctrec[MAXCT*i + ref] stores a histogram of acquisitions of feature i at continuous time reference ref
// times[t] stores the continuous time at which feature t is acquired in the first simulated run
// betas[t] stores the exit propensity after feature t is acquired in the first simulated run
// route[t] is the feature changed at step t
void GetRoutes(int *matrix, int len, int ntarg, double *ntrans, int *rec, double *mean, double *ctrec, double *times, double *betas, int *route, double BINSCALE, int model)
{
  int run, t;
  double time1;
  int state[len];
  double totrate;
  double rate[len];
  double cumsum[len];
  double r;
  int i, j;
  int startt;
  int checker[ntarg];
  int numhits;
  double tgap;
  double continuoustime;

  for(i = 0; i < ntarg; i++)
    checker[i] = 0;

  for(i = 0; i < len; i++)
    mean[i] = 0;
  
  /* loop through NTRAJ simulated trajectories */
  for(run = 0; run < NTRAJ; run++)
    {
      startt = 0; time1 = 0;

      // start at initial state
      for(i = 0; i < len; i++)
	state[i] = 0;

      // track the (continuous) time elapsed
      // (but continuous time is not interpretable unless the posteriors have been produced in the continuous time paradigm)
      continuoustime = 0;

      // loop through feature acquisitions
      for(t = 0; t < len; t++)
	{
	  totrate = 0;
	  // compute the rate for feature i given the current set of features
	  for(i = 0; i < len; i++)
	    {
	      /* ntrans must be the transition matrix. ntrans[i+i*LEN] is the bare rate for i. then ntrans[j*LEN+i] is the modifier for i from j*/
	      if(state[i] == 0)
		{
		  rate[i] = RetrieveEdge(state, i, ntrans, len, model);
		}
	      else // we've already lost this gene
		rate[i] = 0;

	      // roulette wheel calculations as normal
	      cumsum[i] = (i == 0 ? 0 : rate[i-1]+cumsum[i-1]);
	      totrate += rate[i];
	    }

	  // choose a step
	  for(i = 0; i < len; i++)
	    cumsum[i] /= totrate;
	  r = RND;
	  continuoustime += (1./totrate)*log(1./r);

#ifdef VERBOSE
	  for(i = 0; i < len; i++)
	    printf("%.2f ", cumsum[i]);
	  printf("\n");
#endif

	  r = RND;
	  for(i = 0; i < len-1; i++)
	    {
	      if(cumsum[i] < r && cumsum[i+1] > r) { break; }
	    }

#ifdef VERBOSE
	  printf("Rolled %f, chose %i\n", r, i);
#endif

	  // we've chosen feature i, at ordering t, and a timescale continuoustime
	  state[i] = 1;
	  mean[i] += t;

	  // rec[t*len + i] increments if we acquire feature i at ordering t
	  // ctrec[MAXCT*i + ref] increments if we acquire feature i at ct-reference ref
	  // pay attention here! we scale continuous times by BINSCALE (e.g. x100) to produce a reference that allows sensible storage in an integer-referenced histogram, bounded by 0 and MAXCT (element MAXCT-1 stores the number of cases that exceed this)

  	  rec[t*len+i]++;
	  if(continuoustime*BINSCALE < MAXCT)
	    ctrec[MAXCT*i+((int)(BINSCALE*continuoustime))]++;
	  else
	    ctrec[MAXCT*i + MAXCT-1]++;

	  // sample the statistics of the first simulated run. 
	  if(run == 0)
	    {
	      times[t] = continuoustime;
	      betas[t] = totrate;
	      route[t] = i;
	    }

#ifdef VERBOSE
	  for(i = 0; i < len; i++)
	    printf("%i", state[i]);
	  printf(" (%i)\n", t);
#endif

	}
    }

  for(i = 0; i < len; i++)
    mean[i] /= NTRAJ;

}

// construct labels for different features
// for different specific studies this can be adapted to help output
void Label(char *names, int len)
{
  int i, j;
  FILE *fp;
  
  fp = fopen("trait-names.txt", "r");
  if(fp == NULL)
    {
      printf("Didn't find trait-names.txt, using default labels\n");
  for(i = 0; i < len; i++)
    {
      sprintf(&names[i*FLEN], "feature_%i", i);
    }
    }
  else
    {
      i = 0;
      do{
	fgets(&names[i*FLEN], FLEN, fp);
	for(j = 0; j < FLEN; j++)
	  {
	    if(names[i*FLEN+j] == '\n')
	      names[i*FLEN+j] = '\0';
	  }
	i++;
      }while(!feof(fp));
      fclose(fp);
    }
}

void helpandquit(int debug)
{
  printf("Options [defaults]:\n\n--posterior file.txt\tposteriors file [NA]\n--model N\t\tparameter structure (-1 full, 0-4 polynomial degree) [2]\n--seed N\t\trandom seed [0]\n--sims N\t\tsimulations per posterior sample [10]\n--trajs N\t\ttrajectories per simulation [100]\n--burnin N\t\tnumber of samples to skip as burn-in [0]\n--period N\t\tnumber of samples to \"thin\" between sims [0]\n--binscale X\t\tscale for time bins [10]\n--label label\t\tset output file label [OBS FILE AND STATS OF RUN]\n--verbose\t\tverbose file output\n--help\t\t\t[show this message]\n\n");
  exit(0);
}


int main(int argc, char *argv[])
{
  FILE *fp;
  int *matrix;
  int len, ntarg;
  double *trans, *ntrans;
  int t;
  int i, j;
  char ch;
  double lik, nlik;
  int *rec, *order;
  double *drec, *sortdrec, *mean;
  int maxt, allruns;
  int seed = 0;
  char str[200];
  char shotstr[200];
  double tmp;
  int change;
  char names[2000];
  int expt;
  int count;
  double *meanstore, *fmeanstore;
  int specref;
  double *ctrec, ctnorm;
  double *times, *betas;
  int *route;
  FILE *fp1, *fp2, *fp3;
  char fstr[200];
  int tlen;
  int verbose;
  int fref;
  double BINSCALE;
  char postfile[1000];
  int filelabel;
  char labelstr[1000];
  int NVAL;
  int model;
  int burnin, sampleperiod;
  
  // default values
  BINSCALE = 10;
  verbose = 0;
  filelabel = 0;
  seed = 0;
  model = 2;
  burnin = 0;
  sampleperiod = 0;
  sprintf(postfile, "");
  
    printf("\nHyperTraPS(-CT) posterior analysis\n\n");

    // deal with command-line arguments
   for(i = 1; i < argc; i+=2)
    {
      if(strcmp(argv[i], "--posterior\0") == 0) strcpy(postfile, argv[i+1]);
      else if(strcmp(argv[i], "--label\0") == 0) { filelabel = 1; strcpy(labelstr, argv[i+1]); }
      else if(strcmp(argv[i], "--seed\0") == 0) seed = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--model\0") == 0) model = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--sims\0") == 0) NSAMP = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--trajs\0") == 0) NTRAJ = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--burnin\0") == 0) burnin = atoi(argv[i+1]);
      else if(strcmp(argv[i], "--period\0") == 0) sampleperiod = atoi(argv[i+1]);      
      else if(strcmp(argv[i], "--binscale\0") == 0) BINSCALE = atof(argv[i+1]);
      else if(strcmp(argv[i], "--verbose\0") == 0) { verbose = 1; i--; }
      else if(strcmp(argv[i], "--help\0") == 0) helpandquit(0);
    }

   if(strcmp(postfile, "") == 0)
     {
       printf("*** I need at least a file of posterior samples! ***\n\n");
       helpandquit(0);
     }
   if(model == 0)
     {
       printf("*** Posterior analysis isn't meaningful for a zero-parameter model ***\n\n");
       return 0;
     }

  printf("Verbose flag is %i\n", verbose);
  printf("Bin scale is %f\n", BINSCALE);
  
  // open posterior file to assess format
  fp = fopen(postfile, "r");
  if(fp == NULL)
    {
      printf("Couldn't open posterior file %s\n", postfile);
      return 0;
    }
  tlen = 0;
  do{
    ch = fgetc(fp);
    if(ch == ' ') tlen++;
  }while(ch != '\n' && ch != EOF);
  fclose(fp);

  if(ch == EOF)
    {
      printf("Couldn't find appropriate samples in file %s\n", postfile);
      return 0;
    }

  // figure out if posterior file is presented in L*L format; get L if so
  len = 0;
  for(i = 1; i < 200; i++)
    {
      if(model == 1 && tlen == i) { len = i; break; }
      if(model == 2 && tlen == i*i) { len = i; break; }
      if(model == 3 && tlen == i*i*i) { len = i; break; }
      if(model == 4 && tlen == i*i*i*i) { len = i; break; }
      if(model == -1 && tlen == mypow2(i)) { len = i; break; }
    }
  if(len == 0)
    {
      printf("Given model type %i, couldn't determine number of features from %s, which seems to have %i params per sample\n", model, postfile, tlen);
      return 0;
    }

  printf("Based on %s with %i params per model and model %i, there are %i features\n", postfile, tlen, model, len);


  // initialise and allocate a lot of different arrays to compute and store statistics
  srand48(seed);
  allruns  =0;
  ntarg = 0;
  Label(names, len);
  NVAL = nparams(model, len);
  
  matrix = (int*)malloc(sizeof(int)*10000);
  ctrec = (double*)malloc(sizeof(double)*MAXCT*len);
  times = (double*)malloc(sizeof(double)*len);
  betas = (double*)malloc(sizeof(double)*len);
  route = (int*)malloc(sizeof(int)*len);

  trans = (double*)malloc(sizeof(double)*NVAL); /* transition matrix */
  ntrans = (double*)malloc(sizeof(double)*NVAL);
  rec = (int*)malloc(sizeof(int)*len*len); /* stores step ordering, modified by getlikelihood */
  mean = (double*)malloc(sizeof(double)*len);
  meanstore = (double*)malloc(sizeof(double)*len);
  fmeanstore = (double*)malloc(sizeof(double)*len);
  order = (int*)malloc(sizeof(int)*len);
  drec = (double*)malloc(sizeof(double)*len*len);
  sortdrec = (double*)malloc(sizeof(double)*len*len);

  // initialise
  for(i = 0; i < MAXCT*len; i++)
    ctrec[i] = 0;
  ctnorm = 0;

  for(i = 0; i < len*len; i++)
    rec[i] = 0;

  for(i = 0; i < len; i++)
    fmeanstore[i] = 0;

  if(filelabel == 0)
    sprintf(labelstr, "%s", postfile);

  printf("Output label is %s\n", labelstr);
  
  // set up file outputs
  if(verbose)
    {
      sprintf(fstr, "%s-routes.txt", labelstr);
      fp1 = fopen(fstr, "w");
      sprintf(fstr, "%s-betas.txt", labelstr);
      fp2 = fopen(fstr, "w");
      sprintf(fstr, "%s-times.txt", labelstr);
      fp3 = fopen(fstr, "w");
    }
  
      // try to open this file
      fp = fopen(postfile, "r");
      count = 0;
      printf("Working on file %s\n", postfile);

      // loop through posterior samples in this file
      while(!feof(fp))
	{
	  // read in single posterior sample
	  for(i = 0; i < NVAL; i++)
      	    fscanf(fp, "%lf", &ntrans[i]);
	  
	  // this if statement controls which samples get processed
	  // if we want to include burn-in or subsampling, can put it here
	  if(!feof(fp) && count >= burnin && count % (sampleperiod+1) == 0)
	    {
	      // loop through iterations
	      for(j = 0; j< NSAMP; j++)
		{
		  for(i = 0; i < len; i++)
		    meanstore[i] = 0;
		  // simulate behaviour on this posterior and add statistics to counts and histograms
		  GetRoutes(matrix, len, ntarg, ntrans, rec, meanstore, ctrec, times, betas, route, BINSCALE, model);
		  for(i = 0; i < len; i++)
		    fmeanstore[i] += meanstore[i];
		  ctnorm += NTRAJ;
		  allruns++;

		  if(verbose)
		    {
		      for(i = 0; i < len; i++)
 		        fprintf(fp1, "%i ", route[i]);
		      for(i = 0; i < len; i++)
			fprintf(fp2, "%.15f ", betas[i]);
		      for(i = 0; i < len; i++)
			fprintf(fp3, "%.3e ", times[i]);
		      fprintf(fp1, "\n");
		      fprintf(fp2, "\n");
		      fprintf(fp3, "\n");
		    }
		}
	    }
	  count++;
	}
      fclose(fp);
      if(verbose)
	{
	  fclose(fp1);
	  fclose(fp2);
	  fclose(fp3);
	} 


      printf("allruns is %i\n", allruns);

  // output various summaries
  for(i = 0; i < len; i++)
    printf("%i %f\n", i, fmeanstore[i]/allruns);

  // compute mean orderings
  // rec[t*len+i] is prob of obtaining i at time t

  for(i = 0; i < len*len; i++)
    drec[i] = (double)rec[i]/(allruns*NTRAJ);

  for(i = 0; i < len; i++)
    {
      mean[i] = 0;
      order[i] = i;
      for(t = 0; t < len; t++)
	mean[i] += t*drec[t*len+i];
    }

  // simple bubble sort orders features by mean acquisition order
  do{
    change = 0;
    for(i = 0; i < len-1; i++)
      {
	if(mean[i] > mean[i+1])
	  {
	    tmp = mean[i]; mean[i] = mean[i+1]; mean[i+1] = tmp;
	    tmp = order[i]; order[i] = order[i+1]; order[i+1] = tmp;
	    change = 1;
	  }
      }
  }while(change == 1);
  seed--;

  // output the set of summary statistics
  // rec[t*len+i] is prob of obtaining i at time t

  // this produces the heatmap of acquisition probability by feature and order
  // outputs both the original feature ordering and the above mean-sorted references
  sprintf(str, "%s-bubbles.csv", labelstr);
  fp = fopen(str, "w");
  fprintf(fp, "Time,ReorderedIndex,OriginalIndex,Name,Probability\n");
  for(t = 0; t < len; t++)
    {
      for(i = 0; i < len; i++)
	fprintf(fp, "%i,%i,%i,%s,%.15f\n", t, i, order[i], &names[FLEN*order[i]], drec[t*len+order[i]]);
      fprintf(fp, "\n");
    }

  // these appended comments give gnuplot commands for axis labels if required: both in original and mean-sorted orderings
  // commented here, as we're shifting to csv format
  /*  fprintf(fp, "# set xtics (");
  for(i = 0; i < len; i++)
    fprintf(fp, "\"%s\" %i%c", &names[FLEN*order[i]], i, (i == len-1 ? ')' : ','));
  fprintf(fp, "\n");
  fprintf(fp, "# default-order set xtics (");
  for(i = 0; i < len; i++)
    fprintf(fp, "\"%s\" %i%c", &names[FLEN*i], i, (i == len-1 ? ')' : ','));
  fprintf(fp, "\n");

  fprintf(fp, "# (");
  for(i = 0; i < len; i++)
    fprintf(fp, "%i, ", order[i]);
  fprintf(fp, ")\n");

  fprintf(fp, "# ");
  for(i = 0; i < len; i++)
    fprintf(fp, "%s %.4f, ", &names[FLEN*order[i]], mean[i]);
  fprintf(fp, ")\n");

  fclose(fp);*/

  // this stores the time histograms associated with acquisition times for each feature
  // remember here that we've scaled by BINSCALE to store in an integer-referenced array (see GetRoutes())
  sprintf(str, "%s-timehists.csv", labelstr);
  fp = fopen(str, "w");
  fprintf(fp, "OriginalIndex,Time,Probability\n");
  for(i = 0; i < len; i++)
    {
      tmp = 0;
      for(j = 0; j < MAXCT; j++)
	{
	  fprintf(fp, "%i,%f,%.6f\n", i, j/BINSCALE, ctrec[MAXCT*i+j]/ctnorm);
	  tmp += ctrec[MAXCT*i+j]*j;
	}
      printf("%i %.4f\n", i, tmp/ctnorm);
      fprintf(fp, "\n");
    }

  return 0;
}
