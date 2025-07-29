#!/bin/bash

#ja=run_perbase.job
ja=run_perbase_splice_spec.job
rm $ja


indir=../bam_files/
outdir=../perbase_splice_spec
mkdir $outdir

for in_bam in $indir/*GRCh38.unspliced.sorted.bam;
do
	fname=$(basename $in_bam  .GRCh38.unspliced.sorted.bam)
	out_fn=$outdir/$fname.unspliced.out
	echo $in_bam $out_fn >> $ja
done


for in_bam in $indir/*GRCh38.spliced.sorted.bam;
do
        fname=$(basename $in_bam  .GRCh38.spliced.sorted.bam)
        out_fn=$outdir/$fname.spliced.out
	echo $in_bam $out_fn >> $ja
done
