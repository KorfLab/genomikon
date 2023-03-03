use strict;
use warnings;

my $sum = "cksum | cut -f1 -d ' '";

my @TEST = (

	{
		name => "default test",
		cli  => "./testing 1 | $sum",
		pass => "565792770",
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
