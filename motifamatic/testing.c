#include "genomikon.h"
#include "dmg.h"
#include "pwm.h"


int main(int argc, char **argv) {

	// two types of testing
	// really long memory tests
	// default functional tests for "make test"

	init_DNTP();
	char *alph = get_alphabet(8);

	if (argc != 2) gkn_exit("usage: %s <int>", argv[0]);
	int test = atoi(argv[1]);

	if (test == 1) {
		// test 1: convert to/from string
		char *s = "ttNGATYTG";
		gkn_pwm m = str2pwm(s, "test-pwm");
		gkn_pwm_write(m, stdout);
	} else if (test == 2) {
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
	}

	// very long memory tests...



	return 0;
}
