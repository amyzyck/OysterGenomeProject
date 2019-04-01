#!/usr/bin/perl
use File::Slurp; 
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

######################################################
#
# File    : extractCoords.pl
# History : 10/16/2018 -- Created by Kevin Freeman (KF)
#
########################################################
#
# This is a short script that uses a series of perl 
# regexes to extract and format the latitude and
# longitude coordinates from a .kml file (output when
# you export a map from google my maps)
#
#######################################################


########### Command line options ######################################################

my $kml;
my $outFile = "borders.outer";

my $usage = "\nUsage: $0 [options]\n 
Options: 

     -kml               KML file to extract coordinates from
                        
     -out               Name of the outfile that the coordinates will be printed to (default: borders.outer)
                        
     -help              Show this message

";

GetOptions(
   'out=s'          => \$outFile,
   'kml=s'          => \$kml,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($kml) {
    die "\nNo file to parse\n$usage";
};


######################################################################################

my $fileString = read_file($kml);

$fileString =~ /<coordinates>(.*?)<\/coordinates>/s; 
my $coords  = $1;
$coords     =~ s/^\s*(.*?),0\s*$/$1/gm; 
$coords     =~s/,/ /gs; 

my $outFh;
unless(open ($outFh, ">", $outFile)) {
	die "Can't open $outFile for writing, $!";
}

print $outFh $coords; 

