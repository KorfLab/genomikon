#include "genomikon.h"
#include "dmg.h"

static char *usage = "usage: dmgtest <number of tests>";

void dump_pwm(gkn_pwm pwm) {
	printf("%s\n", pwm->name);
	for (int i = 0; i < pwm->size; i++) {
		printf("%d", i);
		for (int j = 0; j < 4; j++) {
			printf(" %.3f", pwm->score[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, char **argv) {
	int n = 0;
	int sigs[7] = {4, 5, 8, 9, 15, 19, 25};
	double P1 = 0.970;
	double P2 = 0.480;
	double P3 = 0.333;
	double p1 = 0.900;
	double p2 = 0.450;
	double p3 = 0.300;

	if (argc != 2) {
		fprintf(stderr, "%s\n", usage);
		exit(1);
	}
	
	n = atoi(argv[1]);
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 7; j++) {
			int sig = sigs[j];
			dmgen dmg = dmgen_new(sig, P1, P2, P3, p1, p2, p3);
			for (int w = 1; w < 4; w++) {
				int limit = (int)pow(sig, w);
				for (int x = 0; x < limit; x++) {
					gkn_pwm pwm = num2pwm(dmg, x, w);
					dump_pwm(pwm);
					gkn_pwm_free(pwm);
				}
			}
			dmgen_free(dmg);
		}
	}

	return 0;
}
