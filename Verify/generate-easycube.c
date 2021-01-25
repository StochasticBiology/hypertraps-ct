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

void outputstr(FILE *fp, int ref)
{
  switch(ref)
    {
    case 0: fprintf(fp, "0 0 0\n"); break;
    case 1: fprintf(fp, "0 0 1\n"); break;
    case 2: fprintf(fp, "0 1 0\n"); break;
    case 3: fprintf(fp, "0 1 1\n"); break;
    case 4: fprintf(fp, "1 0 0\n"); break;
    case 5: fprintf(fp, "1 0 1\n"); break;
    case 6: fprintf(fp, "1 1 0\n"); break;
    case 7: fprintf(fp, "1 1 1\n"); break;
    }
}

int main(void)
{
  double edges[_N*_N];
  int i, j;
  int state, oldstate;
  int t;
  double tc, tau;
  double sampled[_NH];
  double cumsum[_NH];
  double chance[_NH];
  double r;
  int n;
  int step;
  int tmpt;
  double total;
  double scores[300];
  int run;
  double samp1, samp2, oldsamp;
  FILE *fp, *fptime;

  srand48(45);

  for(i = 0; i < _N; i++)
    {
      for(j = 0; j < _N; j++)
	{
	  // edge going from i to j
	  edges[i*_N+j] = 0;
	}
    }
  edges[7*_N+6] = 0.1;
  edges[6*_N+4] = 1;
  edges[4*_N+0] = 0.1;

  fp = fopen("synth-easycube.txt", "w");
  for(i = 0; i < _N; i++)
    {
      for(j = 0; j < _N; j++)
	fprintf(fp, "%f\n", edges[i*_N+j]);
    }
  fclose(fp);

  // first get probability of being in state 1 at time t by sampling trajectories
  fp = fopen("synth-easycube-data.txt", "w");
  fptime = fopen("synth-easycube-time.txt", "w");
  for(n = 0; n < 1000; n++)
    {
      tc = 0; state = 7;
      samp1 = RND*20;
      samp2 = samp1+RND*10;

      for(;tc < 10000;)
	{
	  total = 0; cumsum[0] = 0;
	  for(j = 0; j < _N; j++)
	    {
	      chance[j] = edges[state*_N+j];
	      cumsum[j] = (j == 0 ? 0 : cumsum[j-1])+chance[j];
	      total += chance[j];
	    }
	  if(total == 0) tau = 10000;
	  else tau = expsample(total);

	  if(samp1 != -1 && tc+tau > samp1)
	    {
	      oldstate = state;
	      oldsamp = samp1;
	      samp1 = -1;
	    }
	  if(samp2 != -1 && tc+tau > samp2)
	    {
	      outputstr(fp, state);
	      outputstr(fp, oldstate);
	      fprintf(fptime, "%f\n", samp2-oldsamp);
	      samp2 = -1;
	    }
	  printf("%i %i %f\n", n, state, tau);
	  if(total == 0) break;
	  r = RND;
	  for(j = 0; cumsum[j]/total < r; j++);
	  state = j;
	  tc += tau;
	}
    }

  fclose(fp);
  fclose(fptime);
    
  return 0;
}
