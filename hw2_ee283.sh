#!/bin/bash

mkdir -p /pub/$USER/EE283/DNAseq/rawdata /pub/$USER/EE283/ATACseq/rawdata /pub/$USER/EE283/RNAseq/rawdata

# DNAseq
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/DNAseq"
DestDir="/pub/$USER/EE283/DNAseq/rawdata"
FILES="$SourceDir/*"
for f in $FILES
do
   ff=$(basename $f)
   echo "Processing $ff file..."
   ln -s $SourceDir/$ff $DestDir/$ff
done

# ATACseq 
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/ATACseq"
DestDir="/pub/$USER/EE283/ATACseq/rawdata"

cat $SourceDir/README.ATACseq.txt | tail -n +2 | head -n -3 > $DestDir/ATACseq.labels.txt

File="$DestDir/ATACseq.labels.txt"
while read p
do
   echo "${p}"
   barcode=$(echo $p | cut -f1 -d" ")
   genotype=$(echo $p | cut -f2 -d" ")
   tissue=$(echo $p | cut -f3 -d" ")
   bioRep=$(echo $p | cut -f4 -d" ")
   READ1=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R1.fq.gz")
   READ2=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R2.fq.gz")
   
   echo READ1 $READ1
   echo READ2 $READ2
   # Create the symlink name plus path
   savename1="${genotype}_${tissue}_${bioRep}_R1.fq.gz"
   savename2="${genotype}_${tissue}_${bioRep}_R2.fq.gz"
   # Now make the symlink
   ln -s "$READ1" "$DestDir/$savename1"
   ln -s "$READ2" "$DestDir/$savename2"
done < $File

# RNAseq
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq"
DestDir="/pub/$USER/EE283/RNAseq/rawdata"

# Create labels file
tail -n +2 $SourceDir/RNAseq384_SampleCoding.txt > $DestDir/RNAseq.labels.txt

File="$DestDir/RNAseq.labels.txt"
while read p
do
   # Extract information from the labels file
   sampnumb=$(echo $p | cut -f1 -d" ")
   lanecolumn=$(echo $p | cut -f3 -d" ")
   i7index=$(echo $p | cut -f4 -d" ")
   laneno="${lanecolumn:1}"
   plexno="${lanecolumn:3}"
   
   echo sampnumb: $sampnumb
   echo lanecolumn: $lanecolumn
   echo i7index: $i7index
   echo laneno: $laneno
   echo plexno: $plexno
   
   # Define the source paths for R1 and R2 fastq files
   READ1=$(find $SourceDir/RNAseq384plex_flowcell01/Project_plex${plexno}/Sample_${sampnumb} -type f -iname "*_R1_*.fastq.gz")
   READ2=$(find $SourceDir/RNAseq384plex_flowcell01/Project_plex${plexno}/Sample_${sampnumb} -type f -iname "*_R2_*.fastq.gz")
   
   # Define the destination file names
   savename1="${sampnumb}_R1.fastq.gz"
   savename2="${sampnumb}_R2.fastq.gz"
   
   # Create symbolic links
   if [[ -n "$READ1" && -n "$READ2" ]]; then
       ln -s $READ1 $DestDir/$savename1
       ln -s $READ2 $DestDir/$savename2
   else
       echo "Files for sample number $sampnumb not found."
   fi
   
done < $File

