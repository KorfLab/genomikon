use strict;
use warnings;

my $TS = "data/ch.9940.fa";
my $TG = "data/ch.9940.gff3";
my $sum = "cksum | cut -f1 -d ' '";

my @TEST = (

	{
		name => "isoformer vanilla",
		cli  => "./isoformer $TS | $sum",
		pass => "3597393693",
	},

	{
		name => "isoformer w/ all models",
		cli  => "./isoformer $TS --dpwm data/donor.pwm --apwm data/acceptor.pwm --emm data/exon.mm --imm data/intron.mm --elen data/exon.len --ilen data/intron.len | $sum",
		pass => "2301551904",
	},

	{
		name => "isoformer gff introns",
		cli  => "./isoformer $TS --introns $TG | $sum",
		pass => "3459895910",
	},

	{
		name => "isoformer +gff +models",
		cli  => "./isoformer $TS --introns $TG --dpwm data/donor.pwm --apwm data/acceptor.pwm --emm data/exon.mm --imm data/intron.mm --elen data/exon.len --ilen data/intron.len | $sum",
		pass => "1608473097",
	},

	{
		name => "isocounter",
		cli  => "./isocounter $TS | $sum",
		pass => "2518916367",
	},

	{
		name => "isocounter +gff",
		cli  => "./isocounter $TS --introns $TG | $sum",
		pass => "3088433971",
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
