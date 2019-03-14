#!/usr/bin/perl

use File::Basename;
use File::Path (mkpath);
use strict;

my @split_files = <split_transcript_reads/*>;

mkpath('salmon_alignment');

for (@split_files) {
	my $file_name = basename($_);
	
	my $command = "salmon quant -q --validateMappings -i ../../salmon_indexes/Homo_sapiens.GRCh38.cdna.all_index/ -l A -r $_ -o salmon_alignment/$file_name";
	
	system($command);
	# print("$command\n");
}
