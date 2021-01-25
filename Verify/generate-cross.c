#include <stdio.h>

int main(void)
{
  FILE *fp, *fp1;
  int i, j, k, l;
  char fstr[100];
  int maxj;
  
  for(i = 0; i <= 2; i++)
    {
      sprintf(fstr, "synth-cross-samples-%i.txt", i);
      fp = fopen(fstr, "w");
      sprintf(fstr, "synth-cross-times-%i.txt", i);
      fp1 = fopen(fstr, "w");

      switch(i)
	{
	case 0: maxj = 1; break;
	case 1: maxj = 4; break;
	case 2: maxj = 16; break;
	}
      for(j = 0; j < maxj; j++)
	{
	  for(k = 0; k < 5; k++)
	    {
	      for(l = 0; l < 5; l++)
		fprintf(fp, "%i ", (l < k));
	      fprintf(fp, "\n");
	      for(l = 0; l < 5; l++)
		fprintf(fp, "%i ", (l-1 < k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "0.1\n");
	    }
	  for(k = 4; k >= 0; k--)
	    {
	      for(l = 0; l < 5; l++)
		fprintf(fp, "%i ", (l > k));
	      fprintf(fp, "\n");
	      for(l = 0; l < 5; l++)
		fprintf(fp, "%i ", (l+1 > k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "0.1\n");
	    }
	}

      fclose(fp);
      fclose(fp1);
    }
  
  return 0;
}
		  
