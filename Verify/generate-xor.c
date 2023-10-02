#include <stdio.h>
#include <stdlib.h>

#define RND drand48()

int main(void)
{
  FILE *fp;
  int state[4];
  int i, j;
  int expt;
  
  fp = fopen("synth-xor.txt", "w");
  for(expt = 0; expt < 100; expt++)
    {
      for(i = 0; i < 4; i++)
	state[i] = 0;
  for(i = 0; i < 3; i++)
    {
            for(j = 0; j < 4; j++)
	      fprintf(fp, "%i ", state[j]);
	    fprintf(fp, "\n");
      if(state[0] == 0 && state[1] == 0)
	{
	  if(RND < 0.5) state[0] = 1;
	  else state[1] = 1;
	}
      else if(state[0]+state[1] == 1 && state[2] == 0)
	{
	  if(RND < 0.5) { state[0] = state[1] = 1; }
	  else state[2] = 1;
	}
      else if(state[0]+state[1] == 2)
	{
	  state[3] = 1;
	}
            for(j = 0; j < 4; j++)
	      fprintf(fp, "%i ", state[j]);
	    fprintf(fp, "\n");
    }
    }
  fclose(fp);
  
  return 0;
}
