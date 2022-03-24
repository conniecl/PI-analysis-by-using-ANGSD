# PI-analysis-by-using-ANGSD
#### This is the main procedure that I used for pi analysis, mainly conducted by using ANGSD.
##### 1. align the trimmed data by using Bowtie2 (skip)
##### 2. calculate the pi by using angsd
```
#!/bin/bash
BAMLIST=$1
nInd=$2
OUTFILE=$(echo ${BAMLIST} | sed 's/\.bamlist//')
angsd -bam ${BAMLIST}  -doMaf 1 -doMajorMinor 1 -uniqueOnly 1 -minMapQ 30 -minQ 20 -minInd ${nInd} -doSaf 1 -anc Zea_mays.AGPv4.dna.toplevel.fa -GL 2 -out ${OUTFILE} -fold 1 -P 20
realSFS ${OUTFILE}.saf.idx -P 20 > ${OUTFILE}.sfs
angsd -bam ${BAMLIST}  -doMaf 1 -doMajorMinor 1 -uniqueOnly 1 -minMapQ 30 -minQ 20 -minInd ${nInd} -doSaf 1 -anc Zea_mays.AGPv4.dna.toplevel.fa -GL 2 -out ${OUTFILE} -fold 1 -doThetas 1 -pest ${OUTFILE}.sfs -P 20
thetaStat do_stat ${OUTFILE}.thetas.idx -nChr ${nInd} -win 5000 -step 5000 -outnames ${OUTFILE}.thetasWindow5kb
```
##### 3. calculate the average pi (1Mb as windows)
```
#!usr/bin/perl -w
open OUT,">pi_avg.plot" or die "$!";
@file=qw/maize_rand dip hue lux nic per mex par/;
print OUT "chr\tpos\tpi\tpop\tsite\n";
foreach $in(@file)
{
    open IN,"<$in.thetasWindow5kb.pestPG" or die "$!";
    readline IN;
    %site=();%pi=();
    while(<IN>)
    {
        chomp;
        @tmp=split("\t",$_);
        if($tmp[1]=~/^\d/)
        {
            $pos=int($tmp[2]/1000000);
            if($tmp[-1]!=0)
            {
                $site{$tmp[1]}{$pos}+=$tmp[-1];
                $pi{$tmp[1]}{$pos}+=$tmp[4];
            }
        }
    }
    close IN;
    foreach $key1(sort{$a<=>$b} keys %site)
    {
        foreach $key2(sort{$a<=>$b} keys %{$site{$key1}})
        {
            $avg=$pi{$key1}{$key2}/$site{$key1}{$key2};
            print OUT "$key1\t$key2\t$avg\t$in\n";
        }
    }
}
close OUT;
```
##### 4. Convet the plot file
```
#!usr/bin/perl -w
@file=qw/maize_rand par mex hue dip per lux nic/;
foreach $i(0..$#file)
{
    $hash{$file[$i]}=$i;
}
open IN,"<pi_avg.plot" or die "$!";
readline IN;
open OUT,">pi_avg.convert.plot" or die "$!";
$chr=1;$pos=60;
while(<IN>)
{
    chomp;
    @tmp=split("\t",$_);
    if($tmp[0]!=$chr)
    {
        if($tmp[0]==1)
        {
            $pos=60;
        }
        else
        {
            $chr=$tmp[0];
            $pos=$last+1; $pos+=10
            #print OUT "$tmp[0]\t$pos\t0\t$tmp[3]";
            #$pos++;
            #print OUT "$tmp[0]\t$pos\t0\t$tmp[3]";
            #$pos++;
        }
    }
    $tmp[1]=$tmp[1]+$pos;
    print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$hash{$tmp[3]}\n";
    $last=$tmp[1];
}
close IN;
close OUT;
```
##### 5. visualization it (as in figs7)
```r
library(ggplot2)
data<-read.table("pi_avg.convert.plot")
p<-ggplot(data)+geom_rect(aes(xmin=V2-1,xmax=V2,ymin=V5+1,ymax=V5+1.8,fill=V3))+scale_fill_gradient2(low="#FFFFCC", mid="#E31A1C",high="#800026",midpoint = (max(data$V3)/2))+xlim(0,max(data$V2))+ylim(0,10)+theme_classic()+theme(axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+coord_polar()
ggsave(p,file="pi_avg.png",height=10,width=10)
```r

 


 

