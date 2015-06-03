#Step D in the Workflow Diagram
#This perl script should be placed in the same directory as the data files to be
#processed. It takes the CTD, YSI, and thermistor files after they have been
#initially processed with the proprietary software that came with the sensors,
#and creates versions that will auto-open in Matlab, by stripping off the 
#headers (which often have a variable and unpredictable number of lines). 
#It also creates a file called moor-timestamps.txt to tell Matlab the input 
#variables of importance, that either came from the stripped off headers
#(sensor names, starting date, starting time) or are found manually (starting
#and ending scan numbers for the "good" data). 

#!/usr/bin/perl
use strict;
use warnings;

my $timestamps = "stations\tstartdates\tstarttimes\tstationname\tscanstarts\scanstops\n";

#open the current working directory ( "." )
opendir(my $dh, '.') || die("Cant open dir: $!");

while(my $filename = readdir($dh)) {
#for every thermistor file in the directory:
 if ($filename =~ /.*asc/) {

    open my $file, '<', "$filename" or die ("cant open file: $!");
    my $text = '';
    while(defined(my $line = <$file>)) {
       #clean off whitespaces and newlines
       chomp($line);
       #append this line to our text buffer variable.
       #separate lines with newlines like they were
       $text .= "$line\n";
  }
  close $file;

  my @split = split /start sample number = 1\n/, $text;

  my $station = substr $filename, 0, -4;
  my $starttime = index $split[0], "start time";
  my $date = substr $split[0], $starttime+14, 11;
  my $time = substr $split[0], $starttime+27, 8;
  $timestamps .= "$station\t$date\t$time\tT\t1\t\n";

#from: www.somacon.com/p114.php
sub ltrim($);
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}

my $csv = '';
my @datalines = split "\n", $split[1];
foreach my $dataline (@datalines) {
my @commas = split ",", $dataline;
if (scalar @commas =~ 3) {
$csv .= ltrim($commas[1]).",".ltrim($commas[2]).",".ltrim($commas[0])."\n"; }
# asc file of TD
else {
$csv .= ltrim($commas[2]).",".ltrim($commas[3]).",".ltrim($commas[0]).",".ltrim($commas[1])."\n";
}
}
  open my $datafile, '>', $station . '-matlab.txt' or die ("cant open data file: $!");

  print $datafile $csv;
  close $datafile;

} #closes if thermistor loop

#for the YSI file in the directory
 if ($filename =~ /.*RPT/) {
  open my $file, '<', "$filename" or die ("cant open file: $!");

  my $filelines = '';
  while (defined(my $line = <$file>)) {
      $filelines .= "$line";
  }
  close $file;
  my @text = split /\n/, $filelines;

  my $station = substr $filename, 0, -4;
  open my $newfile, '>', $station . '-matlab.txt' or die ("cant open data file: $!");

#no need for $begin, just use print $newfile $i if $i =~ m!^\d{2}/\d{2}/\d{2} \d{2}:\d{2}!; or look at Regexp::Common for predefined date formats. <-- works if data line is good as-is, which in YSI case isn't.

   my $ysidate = '';
   my $ysitime = '';
   foreach my $i (@text) {
#      if ($i =~ m!^\d{2}/\d{2}/\d{2} \d{2}:\d{2}!)
      my $begin = substr $i, 0, 14;
      if ($begin =~ /From .*/) {
        $ysidate = substr $i, 5, 8;
        $ysitime = substr $i, 14, 5;
      }
      if ($begin =~ /..\/..\/.. ..:../) {
         my @dataline = split / +|\/|:/, $i;
         foreach my $j (@dataline) {print $newfile "$j,"; }
         print $newfile "\n";
      }
   }
  close $newfile;

  $timestamps .= "$station\t$ysidate\t$ysitime\tYSI\t1\t\n";

} #closes if YSI loop

#for every CTD file in the directory (both moored and cast)
 if ($filename =~ /.*cnv/) {

    #open the current file
    open my $file, '<', "$filename" or die ("cant open file: $!");
    my $text = '';
    #for each line in the file
    while(defined(my $line = <$file>)) {
       #clean off whitespaces and newlines
       chomp($line);
       #append this line to our text buffer variable.
       #separate lines with newlines like they were
       $text .= "$line\n";
  }
  close $file;

  my @split = split /\*END\*\n/, $text;

# example: start_time = Jan 01 0000 00:00:00
  my $station = substr $filename, 0, -4;
  my $starttime = index $split[0], "start_time";
  my $date = substr $split[0], $starttime+13, 11;
  my $time = substr $split[0], $starttime+25, 8;

if ($filename =~ /.*_dv.cnv/) { 
  $timestamps .= "$station\t$date\t$time\tCTD\t1\t\n"; 
  open my $datafile, '>', $station . '-matlab.txt' or die ("cant open data file: $!");
  print $datafile $split[1];
  close $datafile;
}
if ($filename =~ /.*cast.cnv/) {
  $timestamps .= "$station\t$date\t$time\tcast\t1\t\n"; 
  open my $datafile, '>', $station . '-matlab.txt' or die ("cant open data file: $!");
  print $datafile $split[1];
  close $datafile;
}
if ($filename =~ /.*_TD.*cnv/) {
  $timestamps .= "$station\t$date\t$time\tTD\t1\t\n"; 
  open my $datafile, '>', $station . 'cnv-matlab.txt' or die ("cant open data file: $!");
  print $datafile $split[1];
  close $datafile;
}

} #closes if cnv loop

} #closes while loop

open my $headerinfo, '>', 'moor-timestamps.txt' or die ("cant open data file: $!");
print $headerinfo "$timestamps";
close $headerinfo;

