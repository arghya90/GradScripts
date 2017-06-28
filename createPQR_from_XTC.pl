#!/usr/bin/perl

use warnings;
use File::Spec;
use File::Basename;

my $help = "\n\nThis script takes a gromacs TPR file which contains all the vital parameters for the initial simulation structure. From that it extracts the info about the resname, atomtype, index, charge and radius. Using that it creates a hash table which stores these info about all the atoms. Then browsing through the XTC file, it only grabs the updated coordinates for each atom at each frame. For all of this, it runs a GMX DUMP to obtain the data in text format for perl to parse. Finally, for each of the frames, a PQR file is output.It has a default radius array for now which is compatible with what gromacs uses.\n RUN IN A QSUBGROMACS MODE\nALSO REMEMBER THAT THE NAME TAGS NEED TO BE CHANGED FOR EACH KIND PF GROMACS SYSTEM. EXAMPLE: PROTEIN_CHAIN_A IS USED HERE WHICH MIGHT BE SYSTEM1 IN A DIFFERENT CASE. MAKE SURE YOU LOOK AT THE TOPOL FILEFOR THIS REFERENCE. SOME BUGS WERE REMOVED IN THIS VERSION BUT THERE ARE CHANCES THAST BUGS MIGHT EXIST ESPECIALLY WHEN IT COMES TO KNOWING THE NUMBER OF DIFFERENT MOLECULES IN THE SYSTEM WHICH DECIDE THE SCALAR KEYS OF MONAME_HASH\n\n";

my $syntax = "perl SCRIPT.pl  <TPR FILE>  <XTC  FILE>  <PREFIX>";

if ( scalar @ARGV != 3 )
{
    print "$help";
    die "$syntax\n\n";
}

#intrinsic variabl3es
my $molname = "";
my $moltype = "";
my $prot_moltype = "";
my $found_the_molecule = 0;
my $numAtoms = 0;
my $numRecords_read = -1;
my $charge_arr_ref = [];
my $radius_arr_ref = [];
my $resind_arr_ref = [];
my $atomName_arr_ref = [];
my $atomName = "";
my $resind = -1;
my $str = "";
my $atomRadius = "";
my $frame = -1;
my $x = -1; my $y = -1; my $z = -1;
my $pqrf = "";

my $vdwradii_hashref = {
  "H"  =>   "0.12",
  "C"  =>   "0.17",
  "N"  =>   "0.155",
  "O"  =>   "0.152",
  "F"  =>   "0.147",
  "P"  =>   "0.18",
  "S"  =>   "0.18"};

my %molname_hash;

#instrinsic functions
sub execute {
  local $sep = "=" x 100;
  print $sep."\n";
  local $cmd = shift;
  print "COMMAND : $cmd\n\n";
  system($cmd);
  print $sep."\n";
}

#user inputs
my $xtcf = $ARGV[1];
my $tprf = $ARGV[0];
my $prefix = $ARGV[2];


execute("gmx dump -s $tprf > tprf.txt ");
open(TPRTXT, "< tprf.txt") or die "NO TPR TEXT FILE FOUND.\nGMX RUN PROBABLY FAILED\n\n";
while( <TPRTXT> )
{
  chomp($_);
  $_ =~ s/^\s{1,}//;
  $_ =~ s/\s{1,}$//;

  if ( $_ =~ m/moltype\s+=\s+(\S+)\s+(\S+)/ )
  {
    $molname = $2;
    $moltype = $1;
    $molname_hash{$moltype} = $molname;
    print "Detected molecule : $molname with molnum $moltype\n\n";
    next;
  }

  # "Will only chose SYSTEM1 as of now because of all the systems have that as the name for the protein part\n";
  if ( scalar(keys(%molname_hash)) == 3)
  {
    foreach ( keys %molname_hash )
    {
      if ( $molname_hash{$_} eq "\"Protein_chain_A\"")
      {
        $molname = "Protein_chain_A";
        $prot_moltype = $_;
        print "Found Protein_chain_A\n";
      }
    }
    undef %molname_hash;
  }

  #now go and look for the lines which has all the details
  if ( $_ =~ m/moltype\s+\((\S+)\):/ )
  {
    $found_the_molecule = 1 if ( $1 eq $prot_moltype);
    next;
  }

  if ($found_the_molecule)
  {
    if ($_ =~ /atom\s+\((\S+)\):/ )
    {
      print $_."\n";
	  $numAtoms = $1;
      print "Expecting $numAtoms atoms in the trajectory. Will exit if that is not the case.\n\n";
      $numRecords_read = 0;

    } elsif ( $_ =~ m/^atom\[\s{1,}[0-9]*\]=/) {

      $atomCharge = substr($_, 69, 12);
      $resind = substr($_,124,5);
      $resind_arr_ref->[$numRecords_read] = $resind;
      $charge_arr_ref->[$numRecords_read] = $atomCharge;
      $numRecords_read += 1;

      $numRecords_read = 0 if ($numRecords_read == $numAtoms);

  } elsif ( $_ =~ m/atom\[(\S+)\]=\{name=\"(\S+)\"\}/ ) {
      $atomName = substr $2,0,1;
      $atomName_arr_ref->[$numRecords_read] = $2;
      $radius_arr_ref->[$numRecords_read] = $vdwradii_hashref->{$atomName};
      $numRecords_read += 1;

      last if ($numRecords_read == $numAtoms);
    }
  }

}

# for ( my $i = 0; $i < $numAtoms; $i += 1)
# {
#     $resind = $resind_arr_ref->[$i];
#     $atomCharge = $charge_arr_ref->[$i];
#     $atomRadius = $radius_arr_ref->[$i];
#     $atomName = $atomName_arr_ref->[$i];
#
#     $str = sprintf("%d\t%4s\t%5d\t%9.4f\t%.4f\n",$i,$atomName,$resind,$atomCharge,$atomRadius);
#     print $str;
# }

close TPRTXT;
execute("rm tprf.txt");

# DONE WITH TPR FILE RECORD READING. ITS TIME TO READ THE DUMPED VERSION OF XTC FILE. IT BETTER BE THE PROTEIN XTC !!!!

execute("gmx dump -f $xtcf > xtcf.txt");
open( XTCTXT, "< xtcf.txt ") or die "NO XTC TEXT FILE FOUND.\nGMX RUN PROBABLY FAILED\n\n";

$numRecords_read = $numAtoms;

while ( <XTCTXT> )
{
    chomp($_);
    $_ =~ s/^\s{1,}//;
    $_ =~ s/\s{1,}$//;

    if ( $_ =~ m/^$xtcf\s+frame\s+(\S+):/)
    {
        $frame = $1;
        $pqrf = $prefix."_frame".$frame.".pqr";

        if ( $numRecords_read == $numAtoms)
        {
            $numRecords_read = 0;
            print "Reading frame $frame... Writing output to $pqrf\n";
            close PQRF if defined PQRF;
            open ( PQRF , "> $pqrf") or warn "FAILED TO WRITE OUTPUT PQR FOR FRAME $frame !!!\n\n";
        }


    } elsif ( ($_ =~ m/^x\[/) ) {
        $resind = $resind_arr_ref->[$numRecords_read];
        $atomName = $atomName_arr_ref->[$numRecords_read];
        $atomCharge = $charge_arr_ref->[$numRecords_read];
        $atomRadius = $radius_arr_ref->[$numRecords_read];

        $x = sprintf "%9.4f" , substr($_,10,12);
        $y = sprintf "%9.4f" , substr($_,23,13);
        $z = sprintf "%9.4f" , substr($_,37,13);

        # print $_."\n";
        # $str = sprintf("%6s%5d%5s%5d%9.4f%9.4f%9.4f%9.4f%9.4f\n","ATOM",$numRecords_read+1,$atomName,$resind,$x*10,$y*10,$z*10,$atomCharge,$atomRadius*10);
        $str = sprintf("%-6s%-5d %-4s%1s%3s %1s%4d   %8.3f%8.3f%8.3f%8.4f%7.4f\n","ATOM",$numRecords_read+1,$atomName,"A","???","X",$resind,$x*10,$y*10,$z*10,$atomCharge,$atomRadius*10);
        # "%-6s%-5d%-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3%8.3f%8.3f\n","ATOM",serial,atomname,"X",resname,"X",resind,x,y,z,q,r
        print PQRF $str;

        $numRecords_read += 1;
    }


}

close XTCTXT;
close PQRF if defined PQRF;

execute("rm xtcf.txt");
