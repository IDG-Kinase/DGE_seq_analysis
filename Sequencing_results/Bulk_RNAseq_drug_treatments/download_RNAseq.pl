#!/usr/bin/perl -w

use strict;
use Text::CSV::Simple;

my $parser = Text::CSV::Simple->new;
my @data = $parser->read_file("Bulk_RNAseq_files.csv");

for (@data[1..$#data]) {
  my $command = "fastq-dump --stdout -gzip ${$_}[0] > ${$_}[1].fastq.gz";
  system($command);
}

