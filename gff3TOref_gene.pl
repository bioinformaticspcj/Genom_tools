#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my($gff_file,%genes1,$geneID,$help,$tran_ID,%trans_cds,%trans_exon_str,%trans_exon_end,%genes2,%chr_bin,%tag,@infor,%exon_frames,@cds_status);
my $bin= 10000;
GetOptions(
        "bin"=>\$bin,
        "gff|:s"=>\$gff_file,
	"help|:s"=>\$help,
)
&usage(),if($help || ! $gff_file);
open GFF3,"$gff_file" or die "can not open $gff_file \n";
my $tag = 0;
$cds_status[0] = "cmpl";
$cds_status[1] = "cmpl";
while(<GFF3>){
	chomp;
	next,if(/^#/);
	my @line = split/\t/;
	$chr_bin{$line[0]} = $bin;
	$tag{$line[0]} = $tag;	
	push @infor,$_;
}
my ($start,$end);
for(@infor){
	chomp;
	my @line = split/\t/;
	$start = $line[3]+1;
	$end = $line[4]+1;
	if($line[2] =~ /gene/i && $_ =~/ID=(\S+?);/){
		$geneID = $1;
		if($line[4] <= $chr_bin{$line[0]} ){
			$genes1{$line[0]}{$geneID} = "$start\t$end\t$line[6]\t$tag{$line[0]}";
		}
		if($line[3] > $chr_bin{$line[0]}){
			$chr_bin{$line[0]} += 10000;
			$tag{$line[0]} += 1;
			$genes1{$line[0]}{$geneID} = "$start\t$end\t$line[6]\t$tag{$line[0]}";
		}
		if($line[3] < $chr_bin{$line[0]} && $line[4] > $chr_bin{$line[0]}){
                        $genes1{$line[0]}{$geneID} = "$start\t$end\t$line[6]\t$tag{$line[0]}";
                }

				
	}
	if($line[2] =~ /mRNA/i && $_ =~/ID=(\S+?);Parent=$geneID;/){
		$tran_ID = $1;
		$genes2{$geneID}{$tran_ID} = "$start\t$end\t$line[6]";
	}
	if($line[2] =~ /exon/i && $_ =~/ID=\S+?;Parent=$tran_ID/){
                $trans_exon_str{$tran_ID}.="$start,";
                $trans_exon_end{$tran_ID}.="$end,";
                $line[7] = "-1",if($line[7] eq ".");
                $exon_frames{$tran_ID}.="$line[7],";
        }
	if($line[2] =~ /CDS/i && $_ =~/ID=\S+?;Parent=$tran_ID/){ 
		if($line[8] =~ /five_prime_partial=true/){
			$cds_status[0] = "incmpl";
		}elsif($line[8] =~ /three_prime_partial=true/){
			$cds_status[1] = "incmpl";
		}else{
			$cds_status[0] = "cmpl";
			$cds_status[1] = "cmpl";
		}		
		$trans_cds{$tran_ID}.="$start,$end,$cds_status[0],$cds_status[1]";	
	}
}
close GFF3;
#open ANO,"$ARGV[1]" or die "can not open $ARGV[1] \n"
foreach my $chr ( keys %genes1){
	foreach my $G_id(keys %{$genes1{$chr}}){
		#print "$G_id\n";
		foreach my $tr_id(keys %{$genes2{$G_id}}){
			my @genes_inf1 = split/\t/,$genes1{$chr}{$G_id};		
			my @genes_inf2 = split/\t/,$genes2{$G_id}{$tr_id};
			my @trans_inf2 = split/,/,$trans_cds{$tr_id};
			my $num_exon =$exon_frames{$tr_id}=~tr/\,/\,/;
			my $out = join"\t",($genes_inf1[3],$G_id,$chr,$genes_inf1[2],$genes_inf2[0],$genes_inf2[1],$trans_inf2[0],$trans_inf2[-3],$num_exon,$trans_exon_str{$tr_id},$trans_exon_end{$tr_id},0,$tr_id,$trans_inf2[2],$trans_inf2[-1],$exon_frames{$tr_id});
			print "$out\n";
		}
	}
}
sub usage{
	my $text=<<"USG";
Name: $0 
	this script can convert a gff3 file to ref_gen file
	-g the input gff3 file
	-b the bin size in ref_gen file (defaut 10000)
	example perl $0 -b 10000 *.gff3 > ref_gen
USG
print $text;
exit;
}
