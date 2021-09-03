use strict;
use warnings;
# use Switch;
use feature qw(switch);

my $fasta_file = $ARGV[0];
my $div = $ARGV[1];
my $fasta_out = 'out_' . $fasta_file;
my $fh_line;
my @fh_header;
my $fh_item;
my $i = 0;
my $fh_seq = '';
my $headline = '';
my $fh_seq_len = 0;
my $q = 0;
my $r = 0;
my $k = 0;
my $l = 0;
my $st = 0;
my $en = 0;
open(FH, "<", $fasta_file) or die "can't open $fasta_file: $!\n";
open(FHO, ">", $fasta_out) or die "can't open $fasta_out: $!\n";
while(<FH>)
{
   next if /^\s*$/;  # skip blank lines
   s/\R//g;
   
   if (/^>/) { # fasta header line
      $fh_line = $_;
      if ($fh_seq ne '') {
         $fh_seq_len = length($fh_seq);
         $headline .= $fh_seq_len . ",";
		 $q = int($fh_seq_len / $div);
		 $r = $fh_seq_len % $div;
		 $st = 0;
		 $en = 0;
		 for( $k = 0; $k < $div; $k++) {
			if ($k < $r) {
				$en = $st + $q + 1;
				$l = $q + 1;
			} else {
				$en = $st + $q;
				$l = $q;
			}
			print FHO "$headline$st,$en,";
            GenNV(substr($fh_seq, $st, $l));
			$st = $en;
		 }
      }
      $fh_seq = '';
      $i++;
      print "$i\n";
      $headline = "$i,";
      $fh_line =~ s/^>gi\|//;
      $fh_line =~ s/Domain\=//;
      @fh_header = split(/\|/, $fh_line);
      foreach my $fh_item (@fh_header) {
         $headline .=  "\"$fh_item\",";
      }
   } else {
      s/\s+//;  # remove any white space
      $fh_seq .= $_;
   }
}

$fh_seq_len = length($fh_seq);
$headline .= $fh_seq_len . ",";
$q = int($fh_seq_len / $div);
$r = $fh_seq_len % $div;
$st = 0;
$en = 0;
for( $k = 0; $k < $div; $k++) {
    if ($k < $r) {
	    $en = $st + $q + 1;
	    $l = $q + 1;
	} else {
	    $en = $st + $q;
	    $l = $q;
	}
	print FHO "$headline$st,$en,";
	GenNV(substr($fh_seq, $st, $l));
	$st = $en;
}



# print FHO length($fh_seq) . ",";
# GenNV($fh_seq);
close(FH);
close(FHO);
1;

sub GenNV {
   my $AminoAcid = 'ACGT';
   my $AA_item;
   my @AA = split(//, $AminoAcid);
   my $SeqText = shift;
   my @Seq = split(//, uc($SeqText));
   my $len = length($SeqText);
   my %n = ();
   my %d = ();
   my %u = ();
   my $i = 0;
   my $j = 0;
   foreach $AA_item (@AA) {
      $n{$AA_item} = 0;
	   $u{$AA_item} = 0;
	   $d{$AA_item} = 0;
#	   $D2{$AA_item} = 0;
#      my @sk = (0) x $len;
#	   $s{$AA_item} = [ @sk ];
   }
   
   for ($i = 0; $i < $len; $i++) {
      if ($Seq[$i] =~ /[ACGT]/) {
         $n{$Seq[$i]}++;
	      $u{$Seq[$i]} += $i + 1;
	      $d{$Seq[$i]} += ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'R') {
         $n{'G'} += 0.5;
         $u{'G'} += 0.5 * ($i + 1);
         $d{'G'} += 0.5 * ($i + 1) ** 2; 
         $n{'A'} += 0.5;
         $u{'A'} += 0.5 * ($i + 1);
         $d{'A'} += 0.5 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'Y') {
         $n{'T'} += 0.5;
         $u{'T'} += 0.5 * ($i + 1);
         $d{'T'} += 0.5 * ($i + 1) ** 2; 
         $n{'C'} += 0.5;
         $u{'C'} += 0.5 * ($i + 1);
         $d{'C'} += 0.5 * ($i + 1) ** 2;
         next;
      }
      if ($Seq[$i] eq  'M') {
         $n{'A'} += 0.5;
         $u{'A'} += 0.5 * ($i + 1);
         $d{'A'} += 0.5 * ($i + 1) ** 2; 
         $n{'C'} += 0.5;
         $u{'C'} += 0.5 * ($i + 1);
         $d{'C'} += 0.5 * ($i + 1) ** 2; 
      }
      if ($Seq[$i] eq  'K') {
         $n{'G'} += 0.5;
         $u{'G'} += 0.5 * ($i + 1);
         $d{'G'} += 0.5 * ($i + 1) ** 2; 
         $n{'T'} += 0.5;
         $u{'T'} += 0.5 * ($i + 1);
         $d{'T'} += 0.5 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'S') {
         $n{'G'} += 0.5;
         $u{'G'} += 0.5 * ($i + 1);
         $d{'G'} += 0.5 * ($i + 1) ** 2; 
         $n{'C'} += 0.5;
         $u{'C'} += 0.5 * ($i + 1);
         $d{'C'} += 0.5 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'W') {
         $n{'A'} += 0.5;
         $u{'A'} += 0.5 * ($i + 1);
         $d{'A'} += 0.5 * ($i + 1) ** 2; 
         $n{'T'} += 0.5;
         $u{'T'} += 0.5 * ($i + 1);
         $d{'T'} += 0.5 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'H') {
         $n{'A'} += 1.0 / 3;
         $u{'A'} += 1.0 / 3 * ($i + 1);
         $d{'A'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'C'} += 1.0 / 3;
         $u{'C'} += 1.0 / 3 * ($i + 1);
         $d{'C'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'T'} += 1.0 / 3;
         $u{'T'} += 1.0 / 3 * ($i + 1);
         $d{'T'} += 1.0 / 3 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'B') {
         $n{'G'} += 1.0 / 3;
         $u{'G'} += 1.0 / 3 * ($i + 1);
         $d{'G'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'C'} += 1.0 / 3;
         $u{'C'} += 1.0 / 3 * ($i + 1);
         $d{'C'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'T'} += 1.0 / 3;
         $u{'T'} += 1.0 / 3 * ($i + 1);
         $d{'T'} += 1.0 / 3 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq 'V') {
         $n{'A'} += 1.0 / 3;
         $u{'A'} += 1.0 / 3 * ($i + 1);
         $d{'A'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'C'} += 1.0 / 3;
         $u{'C'} += 1.0 / 3 * ($i + 1);
         $d{'C'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'G'} += 1.0 / 3;
         $u{'G'} += 1.0 / 3 * ($i + 1);
         $d{'G'} += 1.0 / 3 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'D') {
         $n{'A'} += 1.0 / 3;
         $u{'A'} += 1.0 / 3 * ($i + 1);
         $d{'A'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'G'} += 1.0 / 3;
         $u{'G'} += 1.0 / 3 * ($i + 1);
         $d{'G'} += 1.0 / 3 * ($i + 1) ** 2; 
         $n{'T'} += 1.0 / 3;
         $u{'T'} += 1.0 / 3 * ($i + 1);
         $d{'T'} += 1.0 / 3 * ($i + 1) ** 2; 
         next;
      }
      if ($Seq[$i] eq  'N') {
         $n{'A'} += 0.25;
         $u{'A'} += 0.25 * ($i + 1);
         $d{'A'} += 0.25 * ($i + 1) ** 2; 
         $n{'C'} += 0.25;
         $u{'C'} += 0.25 * ($i + 1);
         $d{'C'} += 0.25 * ($i + 1) ** 2; 
         $n{'G'} += 0.25;
         $u{'G'} += 0.25 * ($i + 1);
         $d{'G'} += 0.25 * ($i + 1) ** 2; 
         $n{'T'} += 0.25;
         $u{'T'} += 0.25 * ($i + 1);
         $d{'T'} += 0.25 * ($i + 1) ** 2; 
         next;
      }
   }
   
   for ($j = 0; $j < length($AminoAcid); $j++) {
      if ($n{$AA[$j]} > 0) {
         $u{$AA[$j]} = $u{$AA[$j]} / $n{$AA[$j]};
         $d{$AA[$j]} = $d{$AA[$j]} / ($n{$AA[$j]} * $len) - $u{$AA[$j]} ** 2 / $len;
      }
      if ($j == (length($AminoAcid) - 1)) {
         print FHO "$n{$AA[$j]},$u{$AA[$j]},$d{$AA[$j]}\n";
      } else {
         print FHO "$n{$AA[$j]},$u{$AA[$j]},$d{$AA[$j]},";
      }
   }
}
   
