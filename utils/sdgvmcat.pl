#!/usr/bin/perl

use strict;

my $dir=shift @ARGV;

unless ($dir && -e $dir) {
  die ("Directory ? $dir\n");
}


my $rout=$dir."/rout";
mkdir $rout unless -e $rout;
if (`ls -1 $rout`) {
  system ("rm -r $rout/*")==0 or die("error: $?");
}
mkdir $rout."/sums" unless -e $rout."/sums";

my %files;
for (my $i=0; $i<1000; $i++) {
  my %files2=%files;
  my $rdir=$dir."/r".$i;
  print $rdir. "\n" ;
  unless (-e $rdir) {
    print "Dernier repertoire: ".($i-1)."\n";
    last;
  }
  opendir(RDIR,$rdir) or die("Ne peux pas ouvrir $rdir");
  while (my $file=readdir(RDIR)) {
    next unless (-f $rdir."/".$file) && ($file=~/\.dat$/);
    system("cat $rdir/$file >>$rout/$file")==0 or die("error: $?");

    #print $rdir."/".$file."\n";
    if ($i==0) { # securite
      $files{$file}=1;
    } else {
      if (!$files{$file}) {
	die ("Ce fichier n'est pas dans les repertoires precedents $rdir");
      }
      delete $files2{$file};
    }
  }
  if (%files2) {
    die("Fichiers introuvables dans $rdir: ".join(" ",keys %files2));
  }
}
