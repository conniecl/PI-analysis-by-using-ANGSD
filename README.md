# PI-analysis-by-using-ANGSD
#### This is the main procedure that I used for pi analysis, mainly conducted by using ANGSD.
##### 1. align the trimmed data by using Bowtie2 (skip)
##### 2. calculate the pi by using angsd (parviglumis, mexicana, nicaraguensis and maize)
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
@file=qw/maize_rand mex par nic/;
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
    foreach $key1(keys %site)
    {
        foreach $key2(keys %{$site{$key1}})
        {
            $avg=$pi{$key1}{$key2}/$site{$key1}{$key2};
            print OUT "chr$key1\t$key2\t$avg\t$in\t$site{$key1}{$key2}\n";
        }
    }
}
close OUT;
```
##### 4. visualization it (as in figs7)
```
library(ggplot2)
data<-read.table("pi_avg.plot",header=T)
data$chr<-factor(data$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"))
p<-ggplot(data)+geom_line(aes(x=pos,y=pi,color=pop))+facet_wrap(~chr,ncol=1,strip.position = "right")+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(size=16),axis.title = element_text(size=21),legend.position=c(0.8,0.2),strip.text=element_text(size=16),strip.background = element_blank())+labs(x="Physical position(Mb)",y="PI")+scale_color_manual(values=c("#0a8646","#1c77c1","#ff932e","#ff3923"))
ggsave(p,filename = "pi_line.svg",height = 12,width =12)
```
##### 5. The relationship between site and pi (use the same input in step3)
```
q<-ggplot(data)+geom_point(aes(x=site,y=pi,color=pop),alpha=1/2)+scale_color_manual(values=c("#0a8646","#1c77c1","#ff932e","#ff3923"))+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(size=16),axis.title = element_text(size=21),legend.position=c(0.8,0.8),strip.text=element_text(size=16),strip.background = element_blank())+labs(x="Site",y="PI")
ggsave(q,filename = "pi_site.png",height = 6,width =6)
```
 


 

