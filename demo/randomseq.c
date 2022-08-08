#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static char *help = "\
usage: randomseq\n\
  -n <int> number of sequences [100]\n\
  -s <int> size of sequences [1000]\n\
  -a <float> fraction of A [0.25]\n\
  -c <float> fraction of C [0.25]\n\
  -g <float> fraction of G [0.25]\n\
  -t <float> fraction of T [0.25]\n\
";

int main(int argc, char **argv) {
	int num = 100;
	int len = 1000;
	double a = 0.25, c = 0.25, g = 0.25, t = 0.25;

	if (argc == 1) {
		printf("%s", help);
		exit(1);
	}

	int opt;
	while ((opt = getopt(argc, argv, "n:s:a:c:g:t:")) != -1) {
		switch (opt) {
			case 'n': num = atoi(optarg); break;
			case 's': len = atoi(optarg); break;
			case 'a': a   = atof(optarg); break;
			case 'c': c   = atof(optarg); break;
			case 'g': g   = atof(optarg); break;
			case 't': t   = atof(optarg); break;
		}
	}

	double total = a + c + g + t;
	a = a / total;
	c = c / total;
	g = g / total;
	t = t / total;

	for (int i = 0; i < num; i++) {
		printf(">seq-%d\n", i);
		for (int j = 0; j < len; j++) {
			double r = random() / (double)RAND_MAX;
			if      (r < a )        putc('A', stdout);
			else if (r < a + c)     putc('C', stdout);
			else if (r < a + c + g) putc('G', stdout);
			else                    putc('T', stdout);

			if (j % 80 == 0 && j != 0) printf("\n");
		}
		printf("\n");
	}
}
