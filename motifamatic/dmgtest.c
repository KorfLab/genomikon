#include "genomikon.h"
#include "dmg.h"



int main(int argc, char **argv) {

	init_DNTP();

	// test 1: convert to/from string
	char *alph = get_alphabet(8);
	char *s = "ACGT";
	gkn_pwm m = str2pwm(s, "test-pwm");
	gkn_pwm_write(m, stdout);

	// test X: make some discretized PWMs of length 3
	int n = 3;
	for (int aid = 0; aid < 8; aid++) {
		int sig = strlen(alph);
		int limit = (int)pow(sig, n);
		for (int x = 0; x < limit; x++) {
			gkn_pwm pwm = num2pwm(alph, x, n);
			gkn_pwm_write(pwm, stdout);
			gkn_pwm_free(pwm);
		}
	}



	return 0;
}
