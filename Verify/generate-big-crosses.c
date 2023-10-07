#include <stdio.h>

int main(void)
{
  FILE *fp, *fp1;
  int i, j, k, l;
  char fstr[100];
  int maxj;
  int L;
  
  for(L = 10; L <= 100; L += 20)
    {
      sprintf(fstr, "synth-bigcross-%i-easy-samples.txt", L);
      fp = fopen(fstr, "w");
      sprintf(fstr, "synth-bigcross-%i-easy-times.txt", L);
      fp1 = fopen(fstr, "w");

   	  for(k = 0; k < L; k++)
	    {
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l < k));
	      fprintf(fp, "\n");
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l-1 < k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "0.1\n");
	    }
	  for(k = L-1; k >= 0; k--)
	    {
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l > k));
	      fprintf(fp, "\n");
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l+1 > k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "0.1\n");
	    }

      fclose(fp);
      fclose(fp1);

            sprintf(fstr, "synth-bigcross-%i-hard-samples.txt", L);
      fp = fopen(fstr, "w");
      sprintf(fstr, "synth-bigcross-%i-hard-times.txt", L);
      fp1 = fopen(fstr, "w");

         	  for(k = 0; k < L; k++)
	    {
	      for(l = 0; l < L; l++)
		fprintf(fp, "0 ");
	      fprintf(fp, "\n");
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l < k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "%f\n", 0.1*k);
	    }
	  for(k = L-1; k >= 0; k--)
	    {
	      for(l = 0; l < L; l++)
		fprintf(fp, "0 ");
	      fprintf(fp, "\n");
	      for(l = 0; l < L; l++)
		fprintf(fp, "%i ", (l > k));
	      fprintf(fp, "\n");
	      fprintf(fp1, "%f\n", 0.1*(L-1-k));
	    }
    }

  return 0;
}
		  
