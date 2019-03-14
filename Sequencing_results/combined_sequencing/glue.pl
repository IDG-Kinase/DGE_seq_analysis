#!/usr/bin/perl -w

use Data::Dumper;
use File::Basename;

my @gz_1 = <'../190205/split_transcript_reads/*'>;
my @gz_2 = <'../190206/split_transcript_reads/*'>;

die if (scalar(@gz_1) != scalar(@gz_2));

for (0..$#gz_1) {
	die if (! basename($gz_1[$_]) eq basename($gz_2[$_]));

	my $command = "cat $gz_1[$_] $gz_2[$_] > 'split_transcript_reads/" . basename($gz_1[$_]) . "'";

	system($command);
}
