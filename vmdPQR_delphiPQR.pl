#!/usr/bin/perl

use warnings;
use File::Basename;

sub print_help {

	my $help = shift;
	my $ip = shift;
	my $op = shift;

	print "PURPOSE : $help";
	print "\n";
	print "INPUT  : $ip";
	print "OUTPUT : $op";
	print "\n";

	die "Exiting ... \n\n";

}

sub VMD2DelPhi_PQR {
	my $help = "This function takes in a PQR file created by VMD or in general any program that generataes a PQR file such that the lines end with a XYZQR format. It converts the file to a PQR format that is compatible with DelPhi\n\n";
	my $isyntax = "PQR FILE FROM VMD\n\n";
	my $osyntax = "PQR FILE COMPATIBLE WITH DELPHI\n\n";

	print_help($help,$isyntax,$osyntax) if (  $_[0] =~ m/\s{0,}/ );
	

	my $atom_key = "";
	my $x = 0;
	my $y = 0;
	my $z = 0;
	my $q = 0;
	my $r = 0;

	my $outline = "";

	my $vmdPQR = $_[0];
	

	my $delphiPQR = basename $vmdPQR;
	
	$delphiPQR =~ s/\.pqr/_delphi\.pqr/;
	$delphiPQR = dirname($vmdPQR)."/".$delphiPQR;

	open ( vPQR, "< $vmdPQR") or die "No such file exists: $vmdPQR \n\n";
	open ( dPQR ,"> $delphiPQR") or die "Could not write the Delphi PQR file\n\n";

	print dPQR "REMARK Converted from VMD to Delphi compatible vers.\n";
	while ( <vPQR> )
	{
		chomp($_);
		
		if ( $_ =~ m/^ATOM/ or $_ =~ m/HETATM/ )
		{
			@line_contents = split(/\s{1,}/ ,$_);
			$r = $line_contents[-1];
			$q = $line_contents[-2];
			$z = $line_contents[-3];
			$y = $line_contents[-4];
			$x = $line_contents[-5];
			$atom_key = substr $_ , 0, 30;

			$outline = sprintf "%s%8.3f%8.3f%8.3f%7.3f%8.4f\n", $atom_key,$x,$y,$z,$q,$r;
			print dPQR "$outline";

		} else {
			print dPQR "$_\n";
		}

	}

	close vPQR;
	close dPQR;
}

my $inPQR = $ARGV[0];
VMD2DelPhi_PQR($inPQR);

