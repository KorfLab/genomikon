#!/usr/bin/perl
use strict;
use warnings;

my @alphabet = qw(A R N D C Q E G H I L K M F P S T W Y V B Z X *);

my %m;
while (<>) {
	next unless $_ =~ /^[A-Z]/;
	chomp;
	my ($let, @val) = split;
	for (my $i = 0; $i < @alphabet; $i++) {
		$m{$let}{$alphabet[$i]} = $val[$i];
	}
}

my @del = qw(B Z *);
foreach my $a1 (@alphabet) {
	foreach my $a2 (@del) {
		delete $m{$a1}{$a2};
	}
}
delete $m{'B'};
delete $m{'Z'};
delete $m{'*'};

my @A = ('A'..'Y');
print "{\n";
foreach my $a1 (@A) {
	my @val;
	foreach my $a2 (@A) {
		if (not defined $m{$a1}{$a2}) {push @val, 0}
		else {push @val, $m{$a1}{$a2}}
	}
	print "\t{", join(",", @val), "},\n";
}
print "};\n";

