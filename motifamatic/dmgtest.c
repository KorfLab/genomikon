#include "genomikon.h"
#include "dmg.h"

static char *usage = "usage: dmgtest <size of motif> <number of tests>";

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
	int sigs[8] = {4, 5, 9, 11, 15, 19, 25, 29};

	if (argc != 3) {
		fprintf(stderr, "%s\n", usage);
		exit(1);
	}

	int s = atoi(argv[1]);
	int n = atoi(argv[2]);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 8; j++) {
			int sig = sigs[j];
			dmgen dmg = dmgen_new(sig);
			int limit = (int)pow(sig, s);
			for (int x = 0; x < limit; x++) {
				gkn_pwm pwm = num2pwm(dmg, x, s);
				dump_pwm(pwm);
				gkn_pwm_free(pwm);
			}
			dmgen_free(dmg);
		}
	}

	return 0;
}
