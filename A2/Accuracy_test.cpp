#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
	int i, count = 0;
	char *gen1, *gen2;
	FILE *f;	
	gen1 = (char *) malloc(1000001);
	gen2 = (char *) malloc(1000001);
	f = fopen(argv[1], "r");
	fscanf(f, "%s", gen1);
	fclose(f);
	f = fopen(argv[2], "r");
	fscanf(f, "%s", gen2);
	fclose(f);
	for (i = 0; i < 1000000; i++) {
		if (gen1[i] == gen2[i]) {
			count++;
		}
	}
	printf("%i\n", count);
	return(0);
}