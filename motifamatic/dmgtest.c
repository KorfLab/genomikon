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
	int sigs[8] = {4, 5, 9, 11, 15, 19, 25, 31};

	if (argc != 2) {
		fprintf(stderr, "%s\n", usage);
		exit(1);
	}

	n = atoi(argv[1]);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 7; j++) {
			int sig = sigs[j];
			dmgen dmg = dmgen_new(sig);
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
