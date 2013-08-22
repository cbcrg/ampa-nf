#!/usr/bin/perl

##
## Multi FASTA wrapper for the AMPA-BIGTABLE.pl core script
##
## If few words, just read a fasta file and for each sequence we obtain a table containing with 6 columns:
    ## Antimicrobial_stretch_Name: Protein_ID.$ ($=number)
    ## Protein_ID
    ## Gene_ID
    ## Start Position of AMP stretch
    ## Length of the AMP stretch (>12aa)
    ## Propensity value of AMP stretch
##

use strict;
use Bio::SeqIO;


my $input_file;
my $result_fileAMP = "AMP_BIGTABLE.txt"; # Main results with one AMP per line
my $result_fileprotein = "Protein_BIGTABLE.txt"; # Data to do a graphic representing the frequency of number of AMP streches per protein
my $target_cmd = '-in=.single.fa.tmp';

my $cl = join( " ", @ARGV);


##
## Print the usage help 
## 
if (($cl=~/-h/) || ($cl=~/-H/) || ($cl=~/-help/) || ($cl=~/-HELP/) ) {
        # print usage 
        print "Usage: multi-AMPA.pl -in=<input fasta file> [other options]\n";
        print "\n";
        print "Available options:\n";
        print "-in      Fasta sequence input file\n";
        print "-t       Threshold value (default: 7)\n";
        print "-w       Window size value (default: 0.225)\n";
        print "-noplot  Skip plot creation step\n";
        print "-help    This help information\n";
        print "\n";
        print "Please note: use the '=' character the separate option name by its value.\n";
        exit 1;
}



##
## option '-in': input file name
##
if ( ($cl=~/-in=\s*(\S+)/)) { 
        $input_file = $1;
}
else {
        print "No input has been specified. Run the script using the \"-in=<file>\" command line option.\n\n";
        exit 1; 
}

# option '-w': window size
if ( ($cl =~ /-w=\s*(\S+)/)) { 
        $target_cmd .= " -w=$1"; 
}

# option '-t': threshold value
if ( ($cl =~ /-t=\s*(\S+)/)) { 
        $target_cmd .= " -t=$1"; 
}


## 
## Read the complete multi-fasta file 
## 
my $in  = Bio::SeqIO->new (
                   -file => $input_file,
                    -format => 'Fasta'
                    );

##
## 
open(AMP, ">$result_fileAMP");
open(PROT, ">$result_fileprotein");

## 
## Iterate over the complete multi-fasta file 
## for each sequence save it in a file named '.single.fa.tmp'
## invoke the AMPA script using that file name 
## 
my $count=0;
my $this_seq;
while ( $this_seq = $in->next_seq() ) {
        $count = $count+1;
        
        open(TMP, ">.single.fa.tmp") or die('Cannot sequence create temporary file');
        print TMP 
                ">", $this_seq->display_id(), "\n",
                ">", $this_seq->desc(), "\n",
                $this_seq->seq(), "\n";
        close(TMP);
        
        my $opt = "$target_cmd -rf=result.$count.txt -df=data.$count.txt " ;

        if ( ($cl =~ /-noplot/)) { 
                $opt .= " -noplot"; 
        }
        else {
                $opt .= " -gf=plot.$count.png "; 
        }
        
        ## 
        ## Invoke the ampa script and check the returned error code
        ## 
        my $out=`perl ./AMPA-BIGTABLE.pl $opt`;
        if ( $? != 0 ) {
          print "Cannot run the AMPA-BIGTABLE.pl script. Returned error: $!\n";
          print "AMPA output: " . $out ."\n"; 
          exit 2;
        }
        
        
        ## 
        ## Appending the result to the single one 
        ## 
        my $seq_desc;
        my @seq_desc;
        my $gene_id;
        my $gene_id1;
        my $protein_id;
        my $prot_aa_sec;
        my $numAMP;
        my $division;
        my $start_pos;

        my $seq_desc = $this_seq->desc();               # Obtaining Gene_ID of each sequence
        my @seq_desc = split( /\s+/, $seq_desc );
        my $gene_id = scalar$seq_desc[2];
        my @gene_id1 = split( /:/, $gene_id );  # Get Gene_ID
        my $protein_id = $this_seq->display_id();   # Get Prot_ID
    
        
        open(RES, "<result.$count.txt") or die('Cannot open the results temporary file');
        while (my $line = <RES>) {  # Read the output file and extract relevant data
            chomp($line);
            if ($line !~ /^\#/) {
                my @result = split( /\s+/, $line );
                my $prot_aa_sec = $result[15];
                my $start_pos = $result[4];
                my $amp_length = $result[6] - $result[4] + 1;   # End - Start Position +1
                my $propensity = $result[9];
                print  AMP $protein_id.".",$count, "\t", $protein_id, "\t", $gene_id1[1], "\t", $start_pos, "\t", $amp_length, "\t", $propensity, "\n"; 
            }
            elsif ($line =~ /^\#/) {
                my @arr = split( /\s+/, $line );
                my $prot_aa_sec = $arr[10];
                my $numAMP = $arr[4];
                my $division = $numAMP/$prot_aa_sec;
                my $amp_aa_length = $arr[15];
                my $best_amp = $arr[19];
                print PROT $protein_id, "\t", $gene_id1[1], "\t", $prot_aa_sec, "\t", $numAMP, "\t", $amp_aa_length, "\t", $best_amp, "\t", $division, "\n";
                    if ($line =~ /0\sbactericidal/ ) {      # If the line contains 0 bactericidal stretches (No AMP)
                        print  AMP $protein_id.".",$count, "\t", $protein_id, "\t", $gene_id1[1], "\t", 0, "\t", 0, "\t", 0, "\n"; 
                    }
            }

        }

        close RES;
        unlink("result.$count.txt");    # Deletes a list of files. 
        unlink("data.$count.txt");
        unlink(".single.fa.tmp");
        
}
close PROT;
close AMP;
