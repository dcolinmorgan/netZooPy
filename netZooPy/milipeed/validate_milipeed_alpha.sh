#!/bin/bash
###BOTH STEPS IN AUROC CALC btwn PWM + wgbsmeth + CHiP meth
## for all ChIP files associated with Motif which will me intersected with WGBS/Me array
#PARAMETER REGION
# validate_milipeed ChIPdir motifdir WGBSdir valoutdir {
# ChIPdir=$1 # %%% location of ChIP bed files
<<<<<<< HEAD
# ChIPdir='data/MotifPipeline/remap/' ## chr start end gene-target per cell line name
# tfdb=$ChIPdir/meta2IDR.txt # metadata file including ChIP file name and gene
# TFfiles=$(ls $ChIPdir*.bed)


# motifdir=$2
motifdir='data/MotifPipeline/hg38_refseq_100kb_tr_fromhg19/'  ## chr start stop pwm per gene name/
# motiffiles=$(ls $motifdir*.txt)

# WGBSdir=$3
WGBSdir='data/MotifPipeline/ENCODE/wgbsin' ## chr start end tmp Me-value, multiple per cell line
=======
ChIPdir='/udd/redmo/data/MotifPipeline/remap/' ## chr start end gene-target per cell line name
# tfdb=$ChIPdir/meta2IDR.txt # metadata file including ChIP file name and gene
TFfiles=$(ls $ChIPdir*.bed)


# motifdir=$2
motifdir='/udd/redmo/data/MotifPipeline/hg38_refseq_100kb_tr_fromhg19/'  ## chr start stop pwm per gene name/
# motiffiles=$(ls $motifdir*.txt)

# WGBSdir=$3
WGBSdir='/udd/redmo/data/MotifPipeline/ENCODE/wgbsin' ## chr start end tmp Me-value, multiple per cell line
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
WGBSfiles=$(ls $WGBSdir/*.bed)
WGBSmeta=$WGBSdir/cellline_meta.txt


# outdir=$4
<<<<<<< HEAD
outdir='data/MotifPipeline/validate_milipeed/'
rm -r -i -f $outdir 
mkdir $outdir
=======
outdir='/udd/redmo/data/MotifPipeline/val_mili_plusNONcpg/'
# rm -r -i -f $outdir
# mkdir $outdir
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
counter=1 
# for wfile in $WGBSfiles
# do
# echo $wfile
<<<<<<< HEAD
# for tfile in $TFfiles # cell line
# do
tfile='data/MotifPipeline/remap/A-549.bed'
=======
for tfile in $TFfiles # cell line
do
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6

TFs=$(eval "cat "$tfile" | cut -f4 | sort | uniq")

for tf in $TFs # targeted genes in ChIP cell lines
do
tftag="${tf}.txt"  #$(eval "echo "$tfile" | cut -d / -f7")

if [ -e $motifdir/$tftag ]
then
<<<<<<< HEAD
# if (( $counter ==1 ))
# then
awk 'NR > 1 || $4 == $tf' $tfile > tffile.txt ## cell-line and gene specific chip file
=======
if (( $counter >250 ))
then
awk 'NR > 1 || $4 == $tf' $tfile > tffilee.txt ## cell-line and gene specific chip file
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6

gtag=$(eval "echo "$tfile" | cut -d / -f7| cut -d . -f1")

###merge multiple wgbs for same cell line (replicates)
<<<<<<< HEAD
# awk -v pat=$gtag '$2 ~ pat'  $WGBSmeta | cut -f 1 | sed 's/^/wgbsin\//' | sed 's/$/.bed/'> tempTF0.txt
# tempTF1=$(cat tempTF0.txt)
=======
awk -v pat=$gtag '$2 ~ pat'  $WGBSmeta | cut -f 1 | sed 's/^/wgbsin\//' | sed 's/$/.bed/'> tempTF00.txt
# tempTF1=$(cat tempTF00.txt)
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
# cat `echo $tempTF1` | cut -f1,2,3,5| sort -V -k1,1 -k2,2n > mergeall.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools merge -i mergeall.txt" > tfile.txt #wgbs file

###use largest of replicates (faster)
<<<<<<< HEAD
# mergeall=$(find $(< tempTF0.txt) | sort -nr | head -n 1 )
=======
mergeall=$(find $(< tempTF00.txt) | sort -nr | head -n 1 )
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6

# cut -f1,2,3,5 $mergall > wfile.txt

echo "$counter : TF=$tf"
counter=$[$counter +1]
### restrict WGBS to only motif regions, return both PWM and %Me
<<<<<<< HEAD
eval "~/bedtools2/bin/bedtools intersect -wa -wb -a $motifdir$tftag -b $WGBSdir/ENCFF005TID.bed" > temp0a.txt
# cut -f1,2,3,4,17 temp0a.txt >temp0b.txt
# min=$(eval "cut -f4 temp0b.txt | scale=4 | bc ")
# min=$(eval "cut -f4 temp0b.txt | sort -n | head -1")
# min=$(cat temp0b.txt | cut -f4 |sort -n | head -1)
# max=$(cat temp0b.txt | cut -f4 |sort -n | tail -1)
# awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),$5/100}' m="$min" ma="$max" OFS='\t' temp0b.txt > temp0d.txt ##standardize range of PWM and Methyl

# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -v -a $motifdir$tftag -b $mergeall" > temp0c.txt
# min=$(cat temp0c.txt | cut -f4 |sort -n | head -1)
# max=$(cat temp0c.txt | cut -f4 |sort -n | tail -1)
# awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),($4-m)/(ma-m)}' m="$min" ma="$max" OFS='\t' temp0c.txt >> temp0d.txt ##combine motif with and without CpG

# cat temp0c.txt |awk 'BEGIN {FS="\t"}; {print $1 $2 $3 $4 $4}' > temp0d.txt
### restrict that intersection above with hits on WGBS, and if no ChIP peak return zero
eval "~/bedtools2/bin/bedtools intersect -wao -a temp0a.txt -b tffile.txt" > $outdir$gtag$tf ##compare entire motif with new methyaltion weights inserted

rm -i -f -r temp0a.txt
# rm -i -f -r temp0b.txt
# rm -i -f -r temp0c.txt
# rm -i -f -r temp0d.txt
rm -i -f -r tempTF0.txt 
rm -i -f -r tffile.txt
# rm -i -f -r sort0bed.txt
# else
# counter=$[$counter +1]
# break
# fi
=======
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $motifdir$tftag -b $mergeall" > temp0aa.txt
cut -f1,2,3,4,17 temp0aa.txt >temp0bb.txt
# min=$(eval "cut -f4 temp0bb.txt | scale=4 | bc ")
# min=$(eval "cut -f4 temp0bb.txt | sort -n | head -1")
min=$(cat temp0bb.txt | cut -f4 |sort -n | head -1)
max=$(cat temp0bb.txt | cut -f4 |sort -n | tail -1)
awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),1-$5/100,1-$5/100}' m="$min" ma="$max" OFS='\t' temp0bb.txt > temp0dd.txt ##standardize range of PWM and Methyl

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -v -a $motifdir$tftag -b $mergeall" > temp0cc.txt
min=$(cat temp0cc.txt | cut -f4 |sort -n | head -1)
max=$(cat temp0cc.txt | cut -f4 |sort -n | tail -1)
awk '{print $1,$2,$3,$4,$5,($4-m)/(ma-m),($4-m)/(ma-m),1}' m="$min" ma="$max" OFS='\t' temp0cc.txt >> temp0dd.txt ##combine motif with and without CpG

# cat temp0cc.txt |awk 'BEGIN {FS="\t"}; {print $1 $2 $3 $4 $4}' > temp0dd.txt
### restrict that intersection above with hits on WGBS, and if no ChIP peak return zero
eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a temp0dd.txt -b tffilee.txt" > $outdir${gtag}_${tf} ##compare entire motif with new methyaltion weights inserted

rm -i -f -r temp0aa.txt
rm -i -f -r temp0bb.txt
rm -i -f -r temp0cc.txt
rm -i -f -r temp0dd.txt
rm -i -f -r tempTF00.txt 
rm -i -f -r tffilee.txt
# rm -i -f -r sort0bed.txt
else
counter=$[$counter +1]
# break
fi
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
fi
done
done
# }
