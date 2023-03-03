#include "genomikon.h"
#include "dmg.h"
#include "pwm.h"

void private_test(void);
void dmg_memory_test(void);
void dmg_test(void);
void pwm_test(void);

int main(int argc, char **argv) {

	init_DNTP();

	if (argc != 2) gkn_exit("usage: %s <int>", argv[0]);
	int test = atoi(argv[1]);
	
	if      (test == 0) private_test();
	else if (test == 1) dmg_memory_test();
	else if (test == 2) dmg_test();
	//else if (test == 3) pwm_test();
	



	return 0;
}

void private_test(void) {
	// whatever I happen to be working on
	gkn_pwm pwm = str2pwm("ttNGATYTG", "IME-motif");
	gkn_pwm_write(pwm, stdout);
	char *str = pwm2str(pwm);
	printf(">>%s<<\n", str);
	gkn_pwm_free(pwm);
}

void dmg_memory_test(void) {
	int len = 2;
	while (1) {
		for (int i = 0; i < 8; i++) {
			char *alph = get_alphabet(i);
			int n = strlen(alph);
			int limit = (int)pow(n, len);
			for (int x = 0; x < limit; x++) {
				gkn_pwm pwm = num2pwm(alph, x, len);
				gkn_pwm_free(pwm);
			}
			free(alph);
		}
	}
}

void dmg_memory_test2(void) {
	// pwm2str
	// str2pwm

}


void dmg_test(void) {
	int len = 2;
	for (int i = 0; i < 8; i++) {
		char *alph = get_alphabet(i);
		int n = strlen(alph);
		int limit = (int)pow(n, len);
		for (int x = 0; x < limit; x++) {
			gkn_pwm pwm = num2pwm(alph, x, len);
			gkn_pwm_write(pwm, stdout);
			gkn_pwm_free(pwm);
		}
	}
}

void pwm_test(void) {

}


/*

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

*/