#!/usr/bin/perl

use strict;

my $nbprocesses=20;

my $count=`wc coords`;
my $nbchunk=($count-1)/$nbprocesses;
$nbchunk=int($nbchunk)+1 if $nbchunk>int($nbchunk);

print STDERR "nbchunk=$nbchunk\n";

if (`ls -1 runoutput`) {
system ("rm -r runoutput/*");
}

my $dirname="run-$$";

for (my $i=0; $i<$nbprocesses; $i++) {
  # write the batch file
  my $runoutput="runoutput/r$i";
  my $run="r$i";
  my $start=$i*$nbchunk+1;
  my $end=($i+1)*$nbchunk;
  $end=$count if $end>$count;

  my $batchfilename="runoutput/batch-$i.dat";
  open (BATCH, ">$batchfilename") or die ("impossible to write \"$batchfilename\
" filename");
  print BATCH '#
#$ -q short.q
#$ -M m.r.lomas@sheffield.ac.uk
cd /home1/bo/bo1mrl/mark
mkdir '.$runoutput.'
echo "host:" $HOSTNAME
echo "dir:" '.$runoutput.'
sdgvm0 input.dat '.$runoutput.' coords '.$start.' '.$end.'

';

  my $continue=1;
  do {
    my $out=`qsub $batchfilename 2>&1`;
    print STDERR $out;
    if ($out=~/job rejected: Only \d+ jobs are allowed per user/) {
      print STDERR "limit for $batchfilename. Retry in 2 minutes\n";
      sleep(60);
    } else {
      $continue=0;
    }
  } while ($continue);
}
