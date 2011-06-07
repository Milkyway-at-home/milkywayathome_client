/* 
   This program takes on stdin a file containing lines in the format

   [NMPZDU](("0x"[0123456789ABCDEFabcdef]){8}" "){4}"#".* "\n"

   e.g.

   N 0xaffe1234 0x5678cafe 0x3FF921FB 0x54442D18 # RN(arccos(0xaffe12345678cafe))

   Let be s the string obtained by concatenating the second and third word without the leading "0x".
   Let x be the 64 bit integer number written in hexadecimal s.

   The program prints out on stdout all lines where x is congruent to 0 modular n 
   The number n is given in argument to the program.

   The program returns 0 if no line has been filtered out.
   It returns 1 if at least one line has been filtered out.
   It returns 2 if it has not been called with a correct argument.
   
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char *argv[]) { 
  int res;
  char *endptr;
  char roundmode;
  unsigned int xh, xl, yh, yl;
  char comment[1024], line[2048];
  unsigned long long int x, n;
  char *l;
  
  if (argc != 2) {
    fprintf(stderr,"Usage: %s n    where n is an integer written in decimal.\n",argv[0]);
    return 2;
  }

  n = (unsigned long long int) strtol(argv[1],&endptr,10);
  
  if (*endptr != '\0') {
    fprintf(stderr,"Usage: %s n    where n is an integer written in decimal.\n",argv[0]);
    return 2;
  }

  res = 0;

  while (!feof(stdin)) {
    if ((fgets(line,2048,stdin) != NULL) && (line[0] != '\n') && (line[0] != '#')) {
   
      sscanf(line,"%c %x %x %x %x %s",&roundmode,&xh,&xl,&yh,&yl,comment);
      
      if (((roundmode != 'N') && 
	   (roundmode != 'U') && 
	   (roundmode != 'D') && 
	   (roundmode != 'Z') &&
	   (roundmode != 'M') &&
	   (roundmode != 'P')) ||
	  ((comment[0] != '\0') &&
	   (comment[0] != '#'))) {
	l = line;
	while ((*l != '\0') && (*l != '\n')) l++;
	if (*l == '\n') *l = '\0';
	fprintf(stderr,"Syntax error in line starting with \"%s\"\n",line);
      } else {
	
	x = xh;
	x <<= 32;
	x += xl;
	
	if ((x % n) == 0) {
	  printf("%s",line);
	} else {
	  res = 1;
	}
      }
    } 
  }

  return res;
}
