#!/usr/bin/perl

use IO::Zlib;
use Data::Dumper;
use Bio::Lite;
use File::Path qw(make_path);
use Text::CSV::Simple;
use strict;

use DateTime;

my $sttime = DateTime->now();

###############################################################################
# Main Program
###############################################################################

my $parser = Text::CSV::Simple->new;
my @data = $parser->read_file('../../plate_ids.csv');

my %read_to_plate;
for (@data) {
	$read_to_plate{$$_[1]} = $$_[0]
}

############################################################
# Use the Map to Sort the Other Sequences
############################################################

$|++;

my $count = 0;
my $matched_seq_count = 0;
my $non_matched_seq_count = 0;

my $barcode_reads = seqFileIterator("All-Sequences_S1_R1_001.fastq.gz");
my $transcript_reads = seqFileIterator("All-Sequences_S1_R2_001.fastq.gz");

my %output_fh;

make_path('split_transcript_reads');

while (my %barcode_entry = %{$barcode_reads->()}) {
	my %transcript_entry = %{$transcript_reads->()};
	
	$barcode_entry{name} =~ /^(NS.*?) /;
	my $barcode_name = $1;
	$transcript_entry{name} =~ /^(NS.*?) /;
	my $transcript_name = $1;
	
	if ($barcode_name != $transcript_name) {
		print "Non-matching entry name: $barcode_name - $transcript_name\n";
	}

	my $this_barcode = substr $barcode_entry{seq}, 0, 6;
	
	$count++;
	if ($count % 1000000 == 0) {
		print "$matched_seq_count $non_matched_seq_count ", $matched_seq_count/$count, " $count ";
		my $entime = DateTime->now;
		my $elapse = $entime - $sttime;
		print "Elapsed time : ".$elapse->in_units('minutes')."m\n";
	}

	if (defined $read_to_plate{$this_barcode}) {
		$matched_seq_count++;
	} else {
		#no well barcode found for this sequnce, break out
		$non_matched_seq_count++;
		next;
	}
	

	my $plate_ID = $read_to_plate{$this_barcode};

	if (! defined $output_fh{$plate_ID}) {
		my $target_file = "split_transcript_reads/$plate_ID.fastq.gz";
		$output_fh{$plate_ID} = IO::Zlib->new($target_file,"rw") or die $!;
	}

	print {$output_fh{$plate_ID}} "\@";
	print {$output_fh{$plate_ID}} "$transcript_entry{name}\n";
	print {$output_fh{$plate_ID}} "$transcript_entry{seq}\n";
	print {$output_fh{$plate_ID}} "+\n";
	print {$output_fh{$plate_ID}} "$transcript_entry{qual}\n";
}
print "$matched_seq_count $non_matched_seq_count ", $matched_seq_count/$count, "\n";

for my $plate_ID (keys %output_fh) {
	undef $output_fh{$plate_ID};
}
