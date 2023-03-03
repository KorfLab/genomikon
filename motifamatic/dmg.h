/******************************************************************************\
 dmg.h
 Copyright (C) Ian Korf
\******************************************************************************/

#ifndef DMG_H
#define DMG_H

#include "genomikon.h"

void    init_DNTP();
void    set_DNTP(double, double, double, double, double, double);
char *  get_alphabet(int);
char *  num2str(const char*, int, int);
gkn_pwm num2pwm(const char*, int, int);
char*   pwm2str(gkn_pwm);
gkn_pwm str2pwm(const char*, const char*);

#endif
