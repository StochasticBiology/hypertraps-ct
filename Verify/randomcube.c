#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _N 8
#define RND drand48()
#define _NH 1000

double expsample(double lambda)
{
  return -log(RND)/lambda;
}

int main(void)
{
  double edges[_N*_N];
  int i, j;
  int state;
  int t;
  double tc, tau;
  double sampled[_NH];
  double cumsum[_NH];
  double chance[_NH];
  double r;
  double v1_path1, v2_path1, v1_path2, v2_path2, w1_path1, w2_path1, w1_path2, w2_path2;
  double u;
  double prob_path1, prob_path2;
  double beta[_NH];
  double analytic;
  double tmpalpha;
  double alpha[_NH];
  double recbeta[_NH*_NH];
  double v1, v2, w1, w2;
  double prob_path;
  int n;
  int step;
  int tmpt;
  double total;
  double scores[300];
  int run;
  FILE *fp;

  srand48(45);

  fp = fopen("randomcubes.txt", "w");
  // loop through 10 random hypercubes
  for(run = 0; run < 10; run++)
    {

      // construct random hypercube
      for(i = 0; i < _N; i++)
	{
	  for(j = 0; j < _N; j++)
	    {
	      // edge going from i to j
	      edges[i*_N+j] = 0;
	      if((i == 7 && (j == 3 || j == 5 || j == 6)) || (i == 3 && (j == 2 || j == 1)) || (i == 5 && (j == 1 || j == 4)) || (i == 6 && (j == 2 || j == 4)) || (i == 1 && j == 0) || (i == 2 && j == 0) || (i == 4 && j == 0))
		edges[i*_N+j] = 0.1+RND;
	    }
	}
      if(run == 0)
	{
	  edges[7*_N+5] = 1;
	  edges[7*_N+3] = 0.1;
	  edges[5*_N+1] = 0.03;
	  edges[3*_N+1] = 10;
	}

      // first get probability of being in state 1 at time t by sampling trajectories

      for(t = 0; t < 100; t++)
	sampled[t] = 0;

      for(n = 0; n < 100000; n++)
	{
	  tc = 0; state = 7;
	  // loop until we get to the opposite corner
	  for(;state != 0;)
	    {
	      // build propensities
	      total = 0; cumsum[0] = 0;
	      for(j = 0; j < _N; j++)
		{
		  chance[j] = edges[state*_N+j];
		  cumsum[j] = (j == 0 ? 0 : cumsum[j-1])+chance[j];
		  total += chance[j];
		}
	      // sample dt
	      tau = expsample(total);
	      if(state == 1 && tc+tau < 100)
		{
		  for(tmpt = ((int)tc)+1; tmpt < tc+tau; tmpt++)
		    sampled[tmpt]++;
		}
	      // choose move
	      r = RND;
	      for(j = 0; cumsum[j]/total < r; j++);
	      state = j;
	      tc += tau;
	    }
	}
      for(t = 0; t < 100; t++)
	{
	  scores[t] = sampled[t]/100000.;
	}

      // now compute it analytically
      for(i = 0; i < _N; i++)
	{
	  beta[i] = 0;
	  for(j = 0; j < _N; j++)
	    beta[i] += edges[i*_N+j];
	}

      // path1 = 111-011-001 7-3-1
      // path2 = 111-101-001 7-5-1
      // easy to compute the explicit probabilities
      prob_path1 = edges[7*_N+3]/beta[7] * edges[3*_N+1]/beta[3];
      prob_path2 = edges[7*_N+5]/beta[7] * edges[5*_N+1]/beta[5];
      v1_path1 = beta[7]*beta[3] * (1./(beta[3]-beta[7]));
      v2_path1 = beta[7]*beta[3] * (1./(beta[7]-beta[3]));
      v1_path2 = beta[7]*beta[5] * (1./(beta[5]-beta[7])); 
      v2_path2 = beta[7]*beta[5] * (1./(beta[7]-beta[5]));
      w1_path1 = beta[7];
      w2_path1 = beta[3];
      w1_path2 = beta[7];
      w2_path2 = beta[5];
      u = -beta[1];

      for(t = 0; t < 100; t ++)
	{
	  analytic = 0;
	  v1 = v1_path1; v2 = v2_path1; w1 = w1_path1; w2 = w2_path1; prob_path = prob_path1;
	  analytic += exp(u*t)*(prob_path*(v1*(1.-exp(-t*(w1 + u)))/(w1 + u) + v2*(1.-exp(-t*(w2 + u)))/(w2 + u)));
	  v1 = v1_path2; v2 = v2_path2; w1 = w1_path2; w2 = w2_path2; prob_path = prob_path2;
	  analytic += exp(u*t)*(prob_path*(v1*(1.-exp(-t*(w1 + u)))/(w1 + u) + v2*(1.-exp(-t*(w2 + u)))/(w2 + u)));
	  scores[100+t] = analytic;
	}

      // now use HyperTraPS, storing betas in recbeta and limiting the choices available by geometry

      for(n = 0; n < _NH; n++)
	{
	  t = 0; state = 7; step = 0; alpha[n] = 1;
	  for(;state != 1;)
	    {
	      tmpalpha = edges[state*_N+5]/beta[state] + edges[state*_N+3]/beta[state] + edges[state*_N+1]/beta[state];
	      recbeta[10*n+step] = beta[state];
	      alpha[n] *= tmpalpha;

	      total = 0; cumsum[0] = 0;
	      for(j = 0; j < _N; j++)
		{
		  chance[j] = (j == 5 || j == 3 || j == 1)*edges[state*_N+j];
		  cumsum[j] = (j == 0 ? 0 : cumsum[j-1])+chance[j];
		  total += chance[j];
		}
	      r = RND;
	      for(j = 0; cumsum[j]/total < r; j++);
	      state = j; step++;
	    }
	}

      for(t = 0; t < 100; t ++)
	{
	  analytic = 0;
	  for(n = 0; n < _NH; n++)
	    {
	      prob_path = alpha[n]*1./_NH;
	      v1 = recbeta[10*n+0]*recbeta[10*n+1] * (1./(recbeta[10*n+1]-recbeta[10*n+0])); 
	      v2 = recbeta[10*n+0]*recbeta[10*n+1] * (1./(recbeta[10*n+0]-recbeta[10*n+1]));
	      w1 = recbeta[10*n+0];
	      w2 = recbeta[10*n+1];
	      u = -beta[1]; 
	      analytic += exp(u*t)*(prob_path*(v1*(1.-exp(-t*(w1 + u)))/(w1 + u) + v2*(1.-exp(-t*(w2 + u)))/(w2 + u)));
	    }
	  scores[200+t] = analytic;
	}
    
      for(i =0; i < 100; i++)
	fprintf(fp, "%i %i %f %f %f\n", run, i, scores[i], scores[100+i], scores[200+i]);
    }
  fclose(fp);

  return 0;
}
