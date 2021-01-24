// this code gets around an awkward formatting issue between NCBI's Newick trees and various Python libraries
// internal nodes need to be labelled for some purposes, so we just create two letter codes "AA", "AB", etc and apply them as we meet internal nodes
// output to stdout
// there is no error checking or anything here -- use with caution!
// also note we are limited to 26^2 internal nodes here!

#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *fp;
  char ch;
  int a;

  // open file
  fp = fopen(argv[1], "r");
  a = 0;

  // loop through tree
  do{
    ch = fgetc(fp);
    if(ch == EOF) break;
    // if we meet a close-paren, append the next unique two-letter code
    if(ch == ')') { printf(")%c%c", 'A'+(a/26), 'A'+(a%26)); a++; fgetc(fp); }
    else printf("%c", ch);
  }while(!feof(fp));
  printf(";\n");

  return 0;
}
