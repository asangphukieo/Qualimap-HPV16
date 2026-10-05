#!/usr/bin/env perl
# bam_coverage.pl - Calculate coverage statistics from samtools mpileup output
# Usage: samtools mpileup sample.bam | perl bam_coverage.pl
# Output: Base_count  Genome_length  Mean_coverage  Min_coverage  Max_coverage

use strict;
use warnings;

my $total_bases = 0;
my $genome_len  = 0;
my $min_cov     = 1e9;
my $max_cov     = 0;

while (<STDIN>) {
    chomp;
    my @fields = split /\t/;
    my $depth  = $fields[3];
    $total_bases += $depth;
    $genome_len++;
    $min_cov = $depth if $depth < $min_cov;
    $max_cov = $depth if $depth > $max_cov;
}

if ($genome_len > 0) {
    my $mean_cov = sprintf("%.2f", $total_bases / $genome_len);
    print "Base_count\tGenome_length\tMean_coverage\tMin_coverage\tMax_coverage\n";
    print "$total_bases\t$genome_len\t$mean_cov\t$min_cov\t$max_cov\n";
} else {
    print "Base_count\tGenome_length\tMean_coverage\tMin_coverage\tMax_coverage\n";
    print "0\t0\t0\t0\t0\n";
}
