#!/usr/bin/perl

# usage: import_data.pl --target-dir ../../raw_data/STORE/2016/2016-05-13_-_Dirnagl_Ulrich/ --link-dir input
# target-dir: path of data folder which should be imported
# link-dir: path of folder which should contain links to 

use strict;
use warnings;
use Cwd qw(getcwd abs_path);
use Getopt::Long;
use File::Basename;

{
    my $wd = getcwd;
    my $target_dir = '';
    my $link_dir = '';
    GetOptions("target-dir=s" => \$target_dir, "link-dir=s" => \$link_dir);
    

    my $base_dir = basename($target_dir);
    my $abs_dir = abs_path($target_dir);

    system("mkdir -p ${abs_dir}");

    my @fastqs = split(/\n/, `ls ${abs_dir}/*/*/*/*.fastq.gz|xargs realpath`);
    
    print join("\t", qw(seq_batch sample flow_cell lane file)),"\n";

    foreach my $target (@fastqs){
	chomp $target;
	my @F = split(/\//, $target);
	my $seq_batch_ix = '';
	foreach my $ix(0..$#F){
	    $seq_batch_ix = $ix if($F[$ix] =~ m/${base_dir}/);
	}
	my $seq_batch = $F[${seq_batch_ix}];

	my $sample = $F[${seq_batch_ix}+1];
	my $flow_cell = $F[${seq_batch_ix}+2];
	my $lane = $F[${seq_batch_ix}+3];
	my $file = $F[${seq_batch_ix}+4];

	print join("\t", ($seq_batch, $sample, $flow_cell, $lane, $file)),"\n";
	my $nd = $link_dir . '/' . $sample;
	system("mkdir -p $nd");
	chdir $nd;
	system("ln -s $target .");
	chdir $wd;
    }
}
