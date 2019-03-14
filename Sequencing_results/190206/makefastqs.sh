#load bcl2fastq2 first
#bsub -n 8 -R "span[hosts=1]" bcl2fastq -p 8 -r 8 -w 8 --no-lane-splitting --runfolder-dir . --output-dir .
sbatch --time=96:00:00 --mem=64G --wrap "bcl2fastq -p 8 -r 8 -w 8 --no-lane-splitting --minimum-trimmed-read-length 17 --mask-short-adapter-reads  17 --runfolder-dir . --output-dir . "

