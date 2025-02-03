#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw( min max );
use Getopt::Long;

# options

## neet to be bam produced by minimap2 with option -<x sr -secondary=yes>
## or equivallently with setting -p -N so that suboptimals are reported
my $inbam = "";

## reads with more than one aloignments with
## edit distance <= (optimalEditDist+$suboptlimit)
## are calssified as multimapper
my $suboptlimit = 1;

## samtools options
my $sto = " --threads 1 ";

GetOptions(
  "i:s"       => \$inbam,
  "l:i"       => \$suboptlimit,
  "st:s"      => \$sto
);

#variables
my $prevqname = "";
my @record    = ();

### add -h below
### remove head
if($inbam =~ m/.sam$/){
  open B, "cat $inbam | " or die "Error: can't open $inbam: $!\n";
  print STDERR "Warning: sam file found. Please take care that sam file is sorted by name."
}
else{
  open B, "samtools sort $sto -n $inbam | samtools view $sto -h  - | " or die "Error: can't open $inbam: $!\n";
}
while(<B>){
  if(m/^@/){
    print;
    next;
  }
  chomp;

  my @F=split"\t",$_;
  if($F[0] eq $prevqname){
    push @record, $_;
    $prevqname = $F[0];
  } else{
    #print STDERR join("\n", @record)."\n";
    &main(\@record) if(@record);
    @record    = ();
    push @record, $_;
    $prevqname = $F[0];
  }
}
close B;

&main(\@record) if(@record);

####################3
# subroutines

sub main{
  my $rec = shift;
  my %reads = ();
  my %shareseq = ();

  foreach my $alnline (@$rec){
    my @F=split"\t",$alnline;

    # skip read if unmapped
    if($F[1] & 4){
      next;
    }

    # skip read if calssified as supplementary alignment
    if($F[1] & 2048){
      next;
    }


    my $qname   = $F[0];
    my $seqname = $F[2];
    my $pos     = $F[3];
    my $mapq    = $F[4];
    my $mate = 0;

    if($F[1] & 64){
      $mate = 1;
    } elsif($F[1] & 128){
      $mate = 2;
    }

    my $editdist = 99;
    if($alnline=~m/NM:i:(\d+)/){
      $editdist = $1;
    }

    push @{$reads{$qname}->{$mate}->{seqname}},    $seqname;
    push @{$reads{$qname}->{$mate}->{pos}},        $pos;
    push @{$reads{$qname}->{$mate}->{mapq}},       $mapq;
    push @{$reads{$qname}->{$mate}->{editdist}},   $editdist;
    push @{$reads{$qname}->{$mate}->{alnline}},    $alnline;

    if($F[9] ne "*"){
      if( defined($shareseq{$qname}->{$mate}->{seq}) ){
        print STDERR "Warning: shareseq already defined for $qname -> $mate -> (<$shareseq{$qname}->{$mate}->{seq}> ~~ <$F[9]>)"."\n";
      }
      $shareseq{$qname}->{$mate}->{seq} = $F[9];
      $shareseq{$qname}->{$mate}->{phred} = $F[10];
    }
  }
  #print Dumper(%reads);


  # report all reads down to optimalEditDist-$suboptlimit
  # report only reads for which mate and pair is reported
  # ?? strict mode: discard reads alltogether if not mate and pair produce the same seq result
  # set mapq=0 if there is more than 1 read reported
  # set NH:i: to number of reads reported

  foreach my $qname (keys %reads){
    my %mateInfo = ();
    my %minEditDistPerSeqname = ();

    foreach my $mate ( keys %{$reads{$qname}} ){
      my $minEditDist = min(@{$reads{$qname}->{$mate}->{editdist}});
      foreach my $idx (  0..$#{$reads{$qname}->{$mate}->{alnline}} ){
        if ($reads{$qname}->{$mate}->{editdist}->[$idx] <= ($minEditDist+$suboptlimit)){
          ##print "$qname -> $mate -> $reads{$qname}->{$mate}->{seqname}->[$idx] -> $reads{$qname}->{$mate}->{editdist}->[$idx]"."\n";
          if( defined($mateInfo{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate}) && defined($minEditDistPerSeqname{$qname}->{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate}) ){
            if( $minEditDistPerSeqname{$qname}->{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate} > $reads{$qname}->{$mate}->{editdist}->[$idx] ){
              $mateInfo{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate} = $idx;
              $minEditDistPerSeqname{$qname}->{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate} = $reads{$qname}->{$mate}->{editdist}->[$idx];
            }
          } else{
            $mateInfo{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate} = $idx;
            $minEditDistPerSeqname{$qname}->{$reads{$qname}->{$mate}->{seqname}->[$idx]}->{$mate} += $reads{$qname}->{$mate}->{editdist}->[$idx];
          }
        }
      }
    }
    #print Dumper(%minEditDistPerSeqname);

    #print Dumper(%mateInfo);
    #print "#####\n";

    my $multimapstatus = 0;
    foreach my $seqname (keys %mateInfo){
      if( defined($mateInfo{$seqname}->{1}) && defined($mateInfo{$seqname}->{2})){
        $multimapstatus++;
      }
    }
    foreach my $seqname (keys %mateInfo){
      if( defined($mateInfo{$seqname}->{1}) && defined($mateInfo{$seqname}->{2})){
        my $idx1   = $mateInfo{$seqname}->{1};
        my $idx2   = $mateInfo{$seqname}->{2};
        my $rname  = $seqname;
        my $pnext1 = $reads{$qname}->{2}->{pos}->[$idx2];
        my $pnext2 = $reads{$qname}->{1}->{pos}->[$idx1];

        #print "$qname -> 1 -> $reads{$qname}->{1}->{seqname}->[$idx1] -> $reads{$qname}->{1}->{editdist}->[$idx1] -> $multimapstatus"."\n";
        #print $reads{$qname}->{1}->{alnline}->[$idx1]."\n";
        print &alnlinemodification($multimapstatus, $rname, $pnext1, $shareseq{$qname}->{1}->{seq}, $shareseq{$qname}->{1}->{phred}, $reads{$qname}->{1}->{alnline}->[$idx1])."\n";

        #print "$qname -> 2 -> $reads{$qname}->{2}->{seqname}->[$idx2] -> $reads{$qname}->{2}->{editdist}->[$idx2] -> $multimapstatus"."\n";
        #print $reads{$qname}->{2}->{alnline}->[$idx2]."\n";
        print &alnlinemodification($multimapstatus, $rname, $pnext2, $shareseq{$qname}->{2}->{seq}, $shareseq{$qname}->{2}->{phred}, $reads{$qname}->{2}->{alnline}->[$idx2])."\n";
      }
    }

  }
}


sub alnlinemodification{
  my $mm = shift;
  my $rn = shift;
  my $pn = shift;
  my $ns = shift;
  my $nq = shift;
  my $li = shift;

  my @li = split"\t", $li;

  if($mm > 1){
    $li[4] = 0;
  }
  $li[6] = $rn;
  $li[7] = $pn;
  if($mm > 1){
    $li[4] = 0;
  }

  if($li[9] eq "*"){
    $li[9] = $ns;
  } else{
    if($li[9] ne $ns){
      print STDERR "Warning: seq does not match (<$li[9]> ~~ <$ns>)\n";
    }
  }

  if($li[10] eq "*"){
    $li[10] = $nq;
  } else{
    if($li[10] ne $nq){
      print STDERR "Warning: phred does not match\n";
    }
  }


  return(join("\t", @li))
}
