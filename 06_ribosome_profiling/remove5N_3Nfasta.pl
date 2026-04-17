#!/usr/bin/perl
#Appends randomers to read name
#Takes one argument, the file name
#Written by Emanuel Wyler, emanuel.wyler@alumni.ethz.ch
use strict;
my $filename =$ARGV[0];
my $N1;
my $N2;
my $filenameout = $filename;
$filenameout =~ s/(.*)\.([^\.]*)$/$1-noUMI\.fasta/;
print ("\nFilename in:",$filename,"\n");
print ("\nFilename out: ",$filenameout,"\n");
open(FILE,$filename); 
open(FOUT,">$filenameout");
while (my $line1 = <FILE>) {
        my $line2 = <FILE>;
        chomp $line1;
        chomp $line2;
        $N1=substr($line2, 0, 3);
        $N2=substr($line2, length($line2)-5, 5);
        $line2 = substr($line2, 3, length($line2)-8);
        print FOUT "$line1;$N1:$N2\n$line2\n";
}
close FILE;
close FOUT;
