use strict;
use warnings;

my $TS = "data/ch.10010.fa";
my $TG = "data/ch.10010.gff3";
my $sum = "cksum | cut -f1 -d ' '";

my @TEST = (

	{
		name => "isoformer vanilla",
		cli  => "./isoformer $TS | $sum",
		pass => "2442921492",
	},

	{
		name => "isoformer w/ all models",
		cli  => "./isoformer $TS --dpwm data/don.pwm --apwm data/acc.pwm --emm data/exon.mm --imm data/intron.mm --elen data/exon.len --ilen data/intron.len | $sum",
		pass => "2196911805",
	},

	{
		name => "isoformer gff introns",
		cli  => "./isoformer $TS --introns $TG | $sum",
		pass => "3336679489",
	},

	{
		name => "isoformer +gff +models",
		cli  => "./isoformer $TS --introns $TG --dpwm data/donor.pwm --apwm data/acceptor.pwm --emm data/exon.mm --imm data/intron.mm --elen data/exon.len --ilen data/intron.len | $sum",
		pass => "2321678176",
	},

	{
		name => "isocounter",
		cli  => "./isocounter $TS | $sum",
		pass => "1826210924",
	},

	{
		name => "isocounter +gff",
		cli  => "./isocounter $TS --introns $TG | $sum",
		pass => "1287712275",
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
