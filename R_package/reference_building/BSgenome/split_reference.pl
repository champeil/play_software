#!/usr/bin/perl                                                                                                                    
# this script is for spliting the reference fa file to each chromosome fa files
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform
# usage
#  enter the result dir and run `perl split_reference.pl [reference file] [output]`
#  then the result will output to the ./

$f = $ARGV[0]; # get the file name
$out = $ARGV[1]; # get the output dir

open (INFILE, "<$f") or die "Can't open: $f $!";

while (<INFILE>) {
        $line = $_; 
        chomp $line;
        if ($line =~ /\>/) { #if > as head
                close OUTFILE;
                my @F = split(/ /,$line);
                my @G = split(/\>/,$F[0]);
                print("$G[1] started \n")
                $line = $F[0];
                $new_file = $out;
                $new_file .= "/";
                $new_file .= $G[1];
                $new_file .= ".fa";
                open (OUTFILE, ">$new_file") or die "Can't open: $new_file $!";
        }
        print OUTFILE "$line\n";
}
close OUTFILE;
