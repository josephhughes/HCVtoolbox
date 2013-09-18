#!/usr/bin/perl

#  this is a perl script for processing a fastq file in combination with the results
#  from a blastn of the fasta/fastq files against a set of reference sequences
#  The fastq file is checked for the presence of barcodes upstream and downstream
#  from the reference match region
#  The fastq region which matches a reference is then checked to correct the sequence
#  if a difference is found relative to the reference. An insert in the reference
#  requires the deletion of the respective site in the fastq, an insert in the
#  fastq (query) requires the insertion of the respective site in the reference
#  Sites that are different to the reference are checked for quality, if the quality is
#  below the threshold and the clean option is specified then the sequence is moved to a notclean file 
# this script is based on fastx toolkit fastx_barcode_splitter from the hannon lab http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
# heavily modified to meet my needs.

use strict;
use Bio::SearchIO;
#use Bio::Index::Fastq;
#use Bio::AlignIO;
use Bio::SeqIO;
use IO::String;
use Getopt::Long; 
use Carp;

# Forward declarations
sub load_barcode_file ($);
sub create_output_files;
sub print_results;
sub open_and_detect_input_format;
sub read_record;
sub match_sequences ;

# Global flags and arguments, 
# Set by command line argumens
my ($blastoutput,$FastqFile,$RefFile,$QualTable,$qthreshold,$outfile,$minlength,$barcode_file,$newfile_prefix,$newfile_suffix);
my $allowed_mismatches = 1;
my $debug = 0 ;
my $barcodes_length;
my $quiet = 0 ;
my $fastq_format = 1;
my $progress=0;
$| = 1;#autoflush stuff

# Global variables 
# Populated by 'create_output_files'
my %barcode;
my %filenames;
my %files;
my %counts = ( 'unmatched' => 0, 'short' => 0 , 'notclean' => 0 );
my %refmatches;
my %uniqseqs;
# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;
my $input_file_io;
my $result;
my $clean;
################



&GetOptions(
	    'inblast:s'     => \$blastoutput,#the blast results
	    'fastq:s'       => \$FastqFile, #fastq file
	    'ref:s'         => \$RefFile, #reference file
	    'qual:s'        => \$QualTable, #quality table for conversion of scores
	    'qcutoff:s'     => \$qthreshold, #quality threshold cut-off
	    'minlength:s'   => \$minlength, #minimum length of a blast match
	    'out:s'         => \$outfile,#output fasta file
	    "bcfile=s"      => \$barcode_file,
	    "prefix=s"      => \$newfile_prefix,
	    "suffix=s"      => \$newfile_suffix,
	    "clean"         => \$clean, #if option clean specified, then sequences that are different to ref with one low quality base are excluded
	    "mismatches=i"  => \$allowed_mismatches,
           );

open(LOG,">Log\_$FastqFile\_$qthreshold\_$minlength\_$allowed_mismatches.txt")||die "Can't open Log.txt\n";
print LOG "inblast $blastoutput fastq $FastqFile ref $RefFile qual $QualTable qcutoff $qthreshold minlength $minlength ";
print LOG "out $outfile bcfile $barcode_file prefix $newfile_prefix suffix $newfile_suffix ";
print LOG "mismatches $allowed_mismatches\n";

if ($clean){
print "clean option selected, only sequences with quality above $qthreshold will be kept\n";
}

#create a loopup table for conversion of quality scores from illumnia 1.3 to phred score
open (QTABLE,"<$QualTable")||die "Can't open $QualTable\n";
my (%quality);
while (<QTABLE>){
  chomp($_);
  my @line=split(/\t/,$_);
  $quality{$line[3]}=$line[0];
}

print "Parsing the blast report...\n";
my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $blastoutput);
# my %blastreport;                           
# while( my $result = $in->next_result ) {
#   ## $result is a Bio::Search::Result::ResultI compliant object  
#   my $query=$result->query_name;
#   print "$query\n";
#   $blastreport{$query}=$result;
# }
# print "BLASTREPORT \n",$blastreport{HLWC18N01DYU7H},"\n";

print "Loading references into memory...\n";
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $RefFile,);
# loads the references sequences into memory - be careful
# if this is a big file, then this script will use a lot of memory
my $refseq;
my %refs;
my %refmatches;
while( $refseq = $seq_in->next_seq() ) {
    $refs{$refseq->display_id()}=$refseq->seq();
}
my $total = `grep -c "^@" $FastqFile 2>&1`;
print "Number of sequences to process $total\n";
print "Beginning the processing ....\n";

open_and_detect_input_format;
load_barcode_file ( $barcode_file ) ;
create_output_files;
match_sequences ;


sub match_sequences {
while ( read_record ) {
  my $percent=$progress*100/$total;
  my $notclean=0;
  printf STDOUT "Progress %03d%%\r",$percent;
  #print "$seq_name$seq_bases$seq_name2$seq_qualities\n";
  chomp($seq_name);
  chomp($seq_bases);
  chomp($seq_name2);
  chomp($seq_qualities);
  if (length($seq_bases)<=187){
    #too short
	my $file = $files{"short"};
   	write_record($file);
   	$counts{'short'}++;
  }
  elsif (length($seq_bases)>187){
  #my $result=$blastreport{$seq_name};
  my $query=$result->query_name;
  die "Error: inconsistent blastreport $query and fastq file $seq_name\n" unless "\@$query" eq $seq_name;
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    print "\nHit=",        $hit->name, " \n";
    my $refname=$hit->name;
    #print $refs{$refname}, "\n";
    print $seq_bases, "\n";
    my ($align1,$align2,$begseq,$endseq)=&SmithWaterman($refs{$refname},$seq_bases);
    
    #print "$begseq\n$endseq\n";
#get the beginning and end and look for the barcode
    $seq_bases=substr ($seq_bases,length($begseq));
    my @newseq=split(//,$seq_bases);
    print "$seq_qualities\n";
    $seq_qualities=substr ($seq_qualities,length($begseq));
    print "\n$align1\n$align2\n\n";#$align1 is the reference
    #print "$seq_qualities\n";
    my @newqual=split(//,$seq_qualities);
    
        #use the unmatched substrings to look for best barcode match
        #first check for exact match on the forward within the first 10 bases
        #save the no matches to unmatched
        #then on the reverse use the match function allowing for 1 mismatch to get the longest substring
        my $best_barcode_mismatches_count = $barcodes_length;
        my $best_barcode_ident = undef;
        my $best_barcode_score=0;
        my $best_insert_length=$minlength;
        foreach my $ident (keys %barcode) {
	      # match forward and reverse primers, ask as input the expected min length of the amplicon
		  my $fbarcode=$barcode{$ident}{"f"};
		  my $rbarcode=	revcompl($barcode{$ident}{"r"});
          my $barcode_score=length($fbarcode)+length($rbarcode);
		  if ($begseq=~/(^.{0,10})$fbarcode/){
		    my $beg_length=length($1);
		    #print "$ident matched in first $beg_length by $fbarcode\n";
		    #print "$endseq $rbarcode $allowed_mismatches\n";
		    my @match= match($endseq, $rbarcode, $allowed_mismatches); #match the longest with  mismatch allowed
		    my $mm = 0;
		    if ($match[0]){
		      $mm=mismatch_count($match[0],$rbarcode);
		    }
		    $barcode_score=$barcode_score-$mm;
		    if ($match[0]){
		      if ($match[0] eq $rbarcode){
		        #print "$ident matched barcode $rbarcode at $match[1] with $match[0] EXACT $match[1]\n";
		      }else{
		        #print "$ident matched barcode $rbarcode at $match[1] with $match[0]\n";
		      }
		      my $insert_length=$match[1];
              if ($best_barcode_score<=$barcode_score){
                  #print "IMMMMMPROVEMENTmatched in first $beg_length by $fbarcode and $ident\t$barcode_score\n";
                  $best_barcode_ident=$ident;
                  $best_barcode_score=$barcode_score;
		      }
		    }
		  }
		}if ($best_barcode_ident){
		  $refmatches{$refname}{$best_barcode_ident}++;
          my $fbarcode=$barcode{$best_barcode_ident}{"f"};
		  my $rbarcode=revcompl($barcode{$best_barcode_ident}{"r"});
		  #print "BEST MACTH for $query is $best_barcode_ident\n";
   		  my $file = $files{$best_barcode_ident};
   		  $counts{$best_barcode_ident}++;
   		  $refmatches{$best_barcode_ident}{$refname}++;
          #now correct the read according to the blast reference match
          #print "QUERY \n$qstr\nHIT   \n$hstr\n";
          #print "query start: $qstart\thit start: $hstart\n";
          my @qstr=split(//,$align2);
          my @hstr=split(//,$align1);
          my $newseqsite;
          my $cnt=0;
          my $mismatches=0;
          for (my $i=0;$i<scalar(@qstr);$i++){
            $newseqsite=$i+$cnt;
            #print "new sequence site: $qstart $newseqsite $newseq[$newseqsite] instead of $i $qstr[$i] hit $hstr[$i]\n";
            if ($qstr[$i] eq $hstr[$i]){
             #print uc($qstr[$i]);
             #print "exact match";
            }elsif ($qstr[$i] ne $hstr[$i]){
              if ($qstr[$i]=~/\-/){
                #print "\nQuery at $i is $qstr[$i] Newseq at: ".($newseqsite)." is ".uc($newseq[$newseqsite])."\nHit at $i is $hstr[$i]\n";
                #check the quality of the bases on either side of this
                #print "=======\nbefore quality ".$quality{$newqual[$newseqsite-1]}." and after ".$quality{$newqual[$newseqsite+1]}."\n=========\n";
                #$newseq[$newseqsite]=~s/(.*)(\w$)/$1$hstr[$i]$2/;
                $newseq[$newseqsite]=~s/(.*)(\w$)/$1N$2/;
                $newqual[$newseqsite]=~s/(.*)(.$)/$1\!$2/;#! phred score of 0 for these types of inserts
                #print "new sequence site: $qstart $newseqsite $newseq[$newseqsite] instead of $i $qstr[$i]\n";
                $cnt--;
                #we need to INSERT this position from the quality line $readq[3] and $readq[1] to keep everything aligned
              }elsif ($hstr[$i]=~/\-/){
                #we need to delete this position from the quality line $readq[3] and $readq[1] to keep everything aligned
                #print "\n>>>>>\nQuery at $i is $qstr[$i] querystart $qstart $hstart Newseq at: ".($newseqsite)." is ".uc($newseq[$newseqsite])."\nHit at $i is $hstr[$i]\n";
                #print "=======\nbefore quality ".$quality{$newqual[$newseqsite-1]}." at site Q ".$quality{$newqual[$newseqsite]}." and after ".$quality{$newqual[$newseqsite+1]}."\n=========\n";
                $newseq[$newseqsite]="";
                $newqual[$newseqsite]="";
              }elsif ($qstr[$i]=~/\w/ && $hstr[$i]=~/\w/){
                #check whether $qstr quality is high enough then print $qstr, else print $hstr
                #print "MISMATCH $i ".uc($qstr[$i])." ".uc($newseq[$newseqsite])." ".$hstr[$i]."\n";
                #print "QUALITY $i ".$newqual[$newseqsite]." ".$quality{$newqual[$newseqsite]}."\n";
                if ($quality{$newqual[$newseqsite]}<$qthreshold){
                  #$newseq[$newseqsite]=$hstr[$i];
                  #$newseq[$newseqsite]="N";
                  if ($clean){#if the clean option has been specified
                    $notclean=1;
                  }  
                  #print "LOW QUALITY at site $newseqsite Q $quality{$newqual[$newseqsite]}\n";
                }else{
                  #print "HIGH QUALITY at site $newseqsite Q $quality{$newqual[$newseqsite]}\n";
                }
              }
            }
          }
          #get the corrected insert sequence and quality
          $seq_bases=join('', @newseq);
          $seq_qualities=join('', @newqual);
          #print "$align1\n$align2\n";
          #print "TRIMMED CORRECTED \n$seq_bases\n$seq_qualities\n";
          $seq_bases=~s/$endseq//;
          my $endlength=length($endseq);
          #$seq_qualities=~s/.{$endlength}$//;
          $seq_qualities=substr ($seq_qualities,1,length($seq_bases));
          print "TRIMMED CORRECTED \n$seq_bases\n$seq_qualities\n";
          if ($notclean=~/1/){
            my $file = $files{"notclean"};
	   	    write_record($file);
	   	    $counts{'notclean'}++;
	   	    $refmatches{$best_barcode_ident}{$refname}--;
	   	    $counts{$best_barcode_ident}--;
	   	    $refmatches{$refname}{$best_barcode_ident}--;
	   	    #print "$seq_name is not clean\n";
          }else{
     	    write_record($file);
     	  }
     	  $uniqseqs{$best_barcode_ident}{$seq_bases}++;
	   	  
	   	}else{
		  #print "NO BARCODE MATCH";
		  #$seq_name=$readq[0];
		  #$seq_bases=$readq[1];
		  #$seq_name2=$readq[2];
		  #$seq_qualities=$readq[3];
   		  my $file = $files{"unmatched"};
	   	  write_record($file);
	   	  $counts{'unmatched'}++;
	   	}
      }
   }  
}
}


print_results unless $quiet;


#
# Read the barcode file
#
sub load_barcode_file ($) {
	my $filename = shift or croak "Missing barcode file name";

	open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
	while (<BCFILE>) {
		next if m/^#/;
		chomp;
		my ($ident, $fbarcode, $rbarcode) = split ;

		$fbarcode = uc($fbarcode);
        $rbarcode = uc($rbarcode);
		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $fbarcode;
		die "Error: bad barcode value ($fbarcode) at barcode file ($filename) line $.\n"
			unless $fbarcode =~ m/^[AGCT]+$/;
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $rbarcode;
		die "Error: bad barcode value ($rbarcode) at barcode file ($filename) line $.\n"
			unless $rbarcode =~ m/^[AGCT]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n" 
			unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident, $fbarcode) is shorter or equal to maximum number of " .
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length($fbarcode)<=$allowed_mismatches;
		die "Error: badcode($ident, $rbarcode) is shorter or equal to maximum number of " .
		    "mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length($rbarcode)<=$allowed_mismatches;
        
       # my $revcompl=revcompl($rbarcode);
        if ($barcode{$ident}{"f"}){
          print "There appears to be multiple forward tags for the identifier $ident\n";
        }if ($barcode{$ident}{"r"}){
          print "There appears to be multiple reverse tags for the identifier $ident\n";
        }if (!$barcode{$ident}{'r'}){
          $barcode{$ident}{'r'}=$rbarcode;
        }if (!$barcode{$ident}{'f'}){
          $barcode{$ident}{'f'}=$fbarcode;
        }
        
	}
	close BCFILE;

	if ($debug) {
		print STDERR "barcode\tforward\treverse\n";
		foreach my $id (keys %barcode) {
			print STDERR $id,"\t", $barcode{$id}{'f'},"\t",$barcode{$id}{'r'},"\n";
		}
	}


}

sub revcompl {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub match {
   my ($s, $t, $max_x) = @_;

   my $m = my @s = unpack('(a)*', $s);
   my $n = my @t = unpack('(a)*', $t);

   push @s, ('?')x($n-1);

   my $best_x = $max_x + 1;
   my $best_i = 0;

   OUTER:
   for my $i (0..$m-1) {
      my $x = 0;

      for my $j (0..$n-1) {
         ++$x if $s[$i+$j] ne $t[$j];
         next OUTER if $x >= $best_x;
      }

      $best_x = $x;
      $best_i = $i;

      last if !$best_x;
   }   

   if ($best_x > $max_x) {
      return undef;
   } else {
      my $mut_barcode = substr($s, $best_i, $n);
      return ($mut_barcode, $best_i);
   }
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub write_record($)
{
	my $file = shift;

	croak "Bad file handle" unless defined $file;

	print $file $seq_name,"\n";
	print $file $seq_bases,"\n";
	#if using FASTQ format, write two more lines
	print $file $seq_name2,"\n";
    print $file $seq_qualities, "\n";
}

sub create_output_files {
    # %barcode should be a unique list of identifiers
	my @exceptions=qw/unmatched short notclean/;
	
   foreach my $except (@exceptions){
     my $new_filename = $newfile_prefix . $except . $newfile_suffix; 
	 $filenames{$except} = $new_filename;
	 open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
	 $files{$except} = $file ;
    }
    
    foreach my $ident (keys %barcode) {
		my $new_filename = $newfile_prefix . $ident . $newfile_suffix; 
		$filenames{$ident} = $new_filename;
		open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n"; 
		$files{$ident} = $file ;
	}
}

sub print_results
{
	print LOG "\nBarcode\tCount\tLocation\tNbDiffRefs";
	for my $refmatch (keys %refmatches){
	  print LOG "$refmatch\t";
	}
	print LOG "\n";
	my $total = 0 ;
	foreach my $ident (sort keys %counts) {
	    my $nbuniqseqs=keys %{$uniqseqs{$ident}};
	    my $nbdifrefmatches=keys %{$refmatches{$ident}};
		print LOG $ident, "\t", $counts{$ident},"\t",$filenames{$ident},"\t",$nbdifrefmatches,"\t";
	    for my $refmatch (keys %refmatches){
	      if ($refmatches{$refmatch}{$ident}){
	        print LOG "$refmatches{$refmatch}{$ident}\t";
	      }else{
	        print LOG "0\t";
	      }
	    }
	    print LOG "\n";
		$total += $counts{$ident};
	}
	print LOG "total\t",$total,"\n";

}

sub open_and_detect_input_format
{
	$input_file_io  = new IO::Handle;
	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");

	# Get the first characeter, and push it back
	my $first_char = $input_file_io->getc();
	$input_file_io->ungetc(ord $first_char);

	if ($first_char eq '>') {
		# FASTA format
		$fastq_format = 0 ;
		print STDERR "Detected FASTA format\n" if $debug;
	} elsif ($first_char eq '@') {
		# FASTQ format
		$fastq_format = 1;
		print STDERR "Detected FASTQ format\n" if $debug;
	} else {
		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
	}
}

sub read_record
{   $result = $in->next_result;
    $progress++;
	$seq_name = $input_file_io->getline();

	return undef unless defined $seq_name; # End of file?

	$seq_bases = $input_file_io->getline();
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;

	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
	}
	return 1;
}

sub SmithWaterman{
my ($seq1,$seq2)=@_;
#print "IN $seq1 $seq2\n";
# scoring scheme
my $MATCH     =  1; # +1 for letters that match
my $MISMATCH = 0; # -1 for letters that mismatch
my $GAP       = -1; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
     $matrix[0][$j]{score}   = 0;
     $matrix[0][$j]{pointer} = "none";
}
for (my $i = 1; $i <= length($seq2); $i++) {
     $matrix[$i][0]{score}   = 0;
     $matrix[$i][0]{pointer} = "none";
}

# fill
 my $max_i     = 0;
 my $max_j     = 0;
 my $max_score = 0;

 for(my $i = 1; $i <= length($seq2); $i++) {
     for(my $j = 1; $j <= length($seq1); $j++) {
         my ($diagonal_score, $left_score, $up_score);
         
         # calculate match score
         my $letter1 = substr($seq1, $j-1, 1);
         my $letter2 = substr($seq2, $i-1, 1);      
         if ($letter1 eq $letter2) {
             $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
          }
         else {
             $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
          }
         
         # calculate gap scores
         $up_score   = $matrix[$i-1][$j]{score} + $GAP;
         $left_score = $matrix[$i][$j-1]{score} + $GAP;
         
         if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
             $matrix[$i][$j]{score}   = 0;
             $matrix[$i][$j]{pointer} = "none";
             next; # terminate this iteration of the loop
          }
         
         # choose best score
         if ($diagonal_score >= $up_score) {
             if ($diagonal_score >= $left_score) {
                 $matrix[$i][$j]{score}   = $diagonal_score;
                 $matrix[$i][$j]{pointer} = "diagonal";
              }
             else {
                 $matrix[$i][$j]{score}   = $left_score;
                 $matrix[$i][$j]{pointer} = "left";
              }
          } else {
             if ($up_score >= $left_score) {
                 $matrix[$i][$j]{score}   = $up_score;
                 $matrix[$i][$j]{pointer} = "up";
              }
             else {
                 $matrix[$i][$j]{score}   = $left_score;
                 $matrix[$i][$j]{pointer} = "left";
              }
          }
         
       # set maximum score
         if ($matrix[$i][$j]{score} > $max_score) {
             $max_i     = $i;
             $max_j     = $j;
             $max_score = $matrix[$i][$j]{score};
          }
      }
 }

 # trace-back

 my $align1 = "";
 my $align2 = "";
 my $foverhang = "";
 my $roverhang = "";
#print "$max_j\n$max_i\n";
 my $j = $max_j;
 my $i = $max_i;
if ($j<$i){
  ###$align1 .= reverse(substr($seq1, $j, length($seq1)-$j));#uncomment ### if you want full alignment with overhangs
  ###$align2 .= reverse(substr($seq2, $i, length($seq2)-$i));
  $roverhang = substr($seq2, $i, length($seq2)-$i);
#  print "Number of gaps to add " , (length($align2)-length($align1)), "\n";
  my $nbgaps=(length($align2)-length($align1));
  if ($nbgaps>0){
    ###$align1 .= "-" x (length($align2)-length($align1));
  }elsif ($nbgaps<0){
    ###$align2 = "-" x (length($align1)-length($align2));
    ###$align2 .= reverse(substr($seq2, $i, length($seq2)-$i));
  }
}

 while (1) {
     if ($matrix[$i][$j]{pointer} eq "none"){
       #print "Min $i, $j\n";
       if ($j<$i){
         ###$align1 .= "-" x $i;
         ###$align2 .= reverse(substr($seq2, $j, $i));
         $foverhang = substr($seq2, $j, $i);
       }elsif ($i<$j){# this hasn't been checked with data yet
         ###$align1 .= reverse(substr($seq2, $i, $j));
         ###$align2 .= "-" x $j;
       }  
       last;
     }  
     if ($matrix[$i][$j]{pointer} eq "diagonal") {#rather than trim substring, I want to do padding
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= substr($seq2, $i-1, 1);
         $i--; $j--;
      }
     elsif ($matrix[$i][$j]{pointer} eq "left") {
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= "-";
         $j--;
      }
     elsif ($matrix[$i][$j]{pointer} eq "up") {
         $align1 .= "-";
         $align2 .= substr($seq2, $i-1, 1);
         $i--;
      }  
 }

 $align1 = reverse $align1;
 $align2 = reverse $align2;
 return ($align1,$align2,$foverhang,$roverhang);
 }