use strict;
use warnings;

my $TF = "values.txt";
my $sum = "cksum | cut -f1 -d ' '";

my @TEST = (

	{
		name => "presti match",
		cli  => "./presti $TF 1 5 20 | $sum",
		pass => "3219087707",
	},

	{
		name => "presti binary",
		cli  => "./presti $TF 1 2 4 8 16 | $sum",
		pass => "3700284193",
	},

);

my $passed = 0;
foreach my $test (@TEST) {
	print STDERR "$test->{name}: ";
	my $result = `$test->{cli}`;
	chomp($result);
	if ($result eq $test->{pass}) {
		$passed++;
		print STDERR "passed\n";
	} else {
		print STDERR "failed\n";
	}
}

print "passed $passed / ", scalar(@TEST), " tests\n";
