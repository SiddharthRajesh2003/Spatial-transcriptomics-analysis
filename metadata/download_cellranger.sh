#!/bin/bash

#SBATCH -J Spatial
#SBATCH -p gpu
#SBATCH -o spatial_%j.txt
#SBATCH -e spatial_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sidrajes@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=18:00:00
#SBATCH --mem=30GB
#SBATCH -A r00750

base=/N/project/Krolab/Siddharth/Pipelines/spatial
cd $base

module load sra-toolkit
# Download and organize FASTQ files by sample
# Files are renamed to Cell Ranger format
# Naming scheme: sequential
# Uses SRA Toolkit (prefetch + fasterq-dump)

set -euo pipefail

# Check if required tools are installed
command -v prefetch >/dev/null 2>&1 || { echo 'prefetch not found. Install SRA Toolkit.'; exit 1; }
command -v fasterq-dump >/dev/null 2>&1 || { echo 'fasterq-dump not found. Install SRA Toolkit.'; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo 'gzip not found. Install gzip for compression.'; exit 1; }

# Sample Mapping:
# GSM6619234 -> primary_tumor_001
# GSM6619236 -> primary_tumor_002
# GSM6619237 -> primary_tumor_003
# GSM6619235 -> recurrent_tumor_001


# ============================================================
# Sample: primary_tumor_001 (Original: GSM6619234)
# Tumor type: primary tumor, Age: 73 years
# ============================================================

echo 'Processing primary_tumor_001...'
mkdir -p fastq/cellranger/primary_tumor_001
cd fastq/cellranger/primary_tumor_001

# Run 1: SRR21832365
echo '  Downloading SRR21832365...'
prefetch SRR21832365 || { echo 'Failed to prefetch SRR21832365'; exit 1; }
fasterq-dump SRR21832365 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832365'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832365_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832365_3.fastq primary_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832365_4.fastq primary_tumor_001_S1_L001_R2_001.fastq
    rm -f SRR21832365_1.fastq SRR21832365_2.fastq
elif [ -f SRR21832365_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832365_2.fastq primary_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832365_3.fastq primary_tumor_001_S1_L001_R2_001.fastq
    rm -f SRR21832365_1.fastq
elif [ -f SRR21832365_1.fastq ] && [ -f SRR21832365_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832365_1.fastq primary_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832365_2.fastq primary_tumor_001_S1_L001_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832365'
    ls -lh SRR21832365*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_001_S1_L001_R1_001.fastq
gzip primary_tumor_001_S1_L001_R2_001.fastq

rm -rf SRR21832365

# Run 2: SRR21832366
echo '  Downloading SRR21832366...'
prefetch SRR21832366 || { echo 'Failed to prefetch SRR21832366'; exit 1; }
fasterq-dump SRR21832366 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832366'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832366_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832366_3.fastq primary_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832366_4.fastq primary_tumor_001_S1_L002_R2_001.fastq
    rm -f SRR21832366_1.fastq SRR21832366_2.fastq
elif [ -f SRR21832366_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832366_2.fastq primary_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832366_3.fastq primary_tumor_001_S1_L002_R2_001.fastq
    rm -f SRR21832366_1.fastq
elif [ -f SRR21832366_1.fastq ] && [ -f SRR21832366_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832366_1.fastq primary_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832366_2.fastq primary_tumor_001_S1_L002_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832366'
    ls -lh SRR21832366*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_001_S1_L002_R1_001.fastq
gzip primary_tumor_001_S1_L002_R2_001.fastq

rm -rf SRR21832366

# Run 3: SRR21832367
echo '  Downloading SRR21832367...'
prefetch SRR21832367 || { echo 'Failed to prefetch SRR21832367'; exit 1; }
fasterq-dump SRR21832367 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832367'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832367_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832367_3.fastq primary_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832367_4.fastq primary_tumor_001_S1_L003_R2_001.fastq
    rm -f SRR21832367_1.fastq SRR21832367_2.fastq
elif [ -f SRR21832367_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832367_2.fastq primary_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832367_3.fastq primary_tumor_001_S1_L003_R2_001.fastq
    rm -f SRR21832367_1.fastq
elif [ -f SRR21832367_1.fastq ] && [ -f SRR21832367_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832367_1.fastq primary_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832367_2.fastq primary_tumor_001_S1_L003_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832367'
    ls -lh SRR21832367*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_001_S1_L003_R1_001.fastq
gzip primary_tumor_001_S1_L003_R2_001.fastq

rm -rf SRR21832367

# Run 4: SRR21832368
echo '  Downloading SRR21832368...'
prefetch SRR21832368 || { echo 'Failed to prefetch SRR21832368'; exit 1; }
fasterq-dump SRR21832368 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832368'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832368_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832368_3.fastq primary_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832368_4.fastq primary_tumor_001_S1_L004_R2_001.fastq
    rm -f SRR21832368_1.fastq SRR21832368_2.fastq
elif [ -f SRR21832368_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832368_2.fastq primary_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832368_3.fastq primary_tumor_001_S1_L004_R2_001.fastq
    rm -f SRR21832368_1.fastq
elif [ -f SRR21832368_1.fastq ] && [ -f SRR21832368_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832368_1.fastq primary_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832368_2.fastq primary_tumor_001_S1_L004_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832368'
    ls -lh SRR21832368*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_001_S1_L004_R1_001.fastq
gzip primary_tumor_001_S1_L004_R2_001.fastq

rm -rf SRR21832368

cd ../../..
echo 'Completed primary_tumor_001'
echo ''


# ============================================================
# Sample: recurrent_tumor_001 (Original: GSM6619235)
# Tumor type: recurrent tumor, Age: 58 years
# ============================================================

echo 'Processing recurrent_tumor_001...'
mkdir -p fastq/cellranger/recurrent_tumor_001
cd fastq/cellranger/recurrent_tumor_001

# Run 1: SRR21832360
echo '  Downloading SRR21832360...'
prefetch SRR21832360 || { echo 'Failed to prefetch SRR21832360'; exit 1; }
fasterq-dump SRR21832360 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832360'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832360_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832360_3.fastq recurrent_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832360_4.fastq recurrent_tumor_001_S1_L001_R2_001.fastq
    rm -f SRR21832360_1.fastq SRR21832360_2.fastq
elif [ -f SRR21832360_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832360_2.fastq recurrent_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832360_3.fastq recurrent_tumor_001_S1_L001_R2_001.fastq
    rm -f SRR21832360_1.fastq
elif [ -f SRR21832360_1.fastq ] && [ -f SRR21832360_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832360_1.fastq recurrent_tumor_001_S1_L001_R1_001.fastq
    mv SRR21832360_2.fastq recurrent_tumor_001_S1_L001_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832360'
    ls -lh SRR21832360*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip recurrent_tumor_001_S1_L001_R1_001.fastq
gzip recurrent_tumor_001_S1_L001_R2_001.fastq

rm -rf SRR21832360

# Run 2: SRR21832361
echo '  Downloading SRR21832361...'
prefetch SRR21832361 || { echo 'Failed to prefetch SRR21832361'; exit 1; }
fasterq-dump SRR21832361 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832361'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832361_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832361_3.fastq recurrent_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832361_4.fastq recurrent_tumor_001_S1_L002_R2_001.fastq
    rm -f SRR21832361_1.fastq SRR21832361_2.fastq
elif [ -f SRR21832361_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832361_2.fastq recurrent_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832361_3.fastq recurrent_tumor_001_S1_L002_R2_001.fastq
    rm -f SRR21832361_1.fastq
elif [ -f SRR21832361_1.fastq ] && [ -f SRR21832361_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832361_1.fastq recurrent_tumor_001_S1_L002_R1_001.fastq
    mv SRR21832361_2.fastq recurrent_tumor_001_S1_L002_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832361'
    ls -lh SRR21832361*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip recurrent_tumor_001_S1_L002_R1_001.fastq
gzip recurrent_tumor_001_S1_L002_R2_001.fastq

rm -rf SRR21832361

# Run 3: SRR21832362
echo '  Downloading SRR21832362...'
prefetch SRR21832362 || { echo 'Failed to prefetch SRR21832362'; exit 1; }
fasterq-dump SRR21832362 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832362'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832362_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832362_3.fastq recurrent_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832362_4.fastq recurrent_tumor_001_S1_L003_R2_001.fastq
    rm -f SRR21832362_1.fastq SRR21832362_2.fastq
elif [ -f SRR21832362_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832362_2.fastq recurrent_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832362_3.fastq recurrent_tumor_001_S1_L003_R2_001.fastq
    rm -f SRR21832362_1.fastq
elif [ -f SRR21832362_1.fastq ] && [ -f SRR21832362_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832362_1.fastq recurrent_tumor_001_S1_L003_R1_001.fastq
    mv SRR21832362_2.fastq recurrent_tumor_001_S1_L003_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832362'
    ls -lh SRR21832362*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip recurrent_tumor_001_S1_L003_R1_001.fastq
gzip recurrent_tumor_001_S1_L003_R2_001.fastq

rm -rf SRR21832362

# Run 4: SRR21832363
echo '  Downloading SRR21832363...'
prefetch SRR21832363 || { echo 'Failed to prefetch SRR21832363'; exit 1; }
fasterq-dump SRR21832363 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832363'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832363_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832363_3.fastq recurrent_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832363_4.fastq recurrent_tumor_001_S1_L004_R2_001.fastq
    rm -f SRR21832363_1.fastq SRR21832363_2.fastq
elif [ -f SRR21832363_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832363_2.fastq recurrent_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832363_3.fastq recurrent_tumor_001_S1_L004_R2_001.fastq
    rm -f SRR21832363_1.fastq
elif [ -f SRR21832363_1.fastq ] && [ -f SRR21832363_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832363_1.fastq recurrent_tumor_001_S1_L004_R1_001.fastq
    mv SRR21832363_2.fastq recurrent_tumor_001_S1_L004_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832363'
    ls -lh SRR21832363*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip recurrent_tumor_001_S1_L004_R1_001.fastq
gzip recurrent_tumor_001_S1_L004_R2_001.fastq

rm -rf SRR21832363

cd ../../..
echo 'Completed recurrent_tumor_001'
echo ''


# ============================================================
# Sample: primary_tumor_002 (Original: GSM6619236)
# Tumor type: primary tumor, Age: 71 years
# ============================================================

echo 'Processing primary_tumor_002...'
mkdir -p fastq/cellranger/primary_tumor_002
cd fastq/cellranger/primary_tumor_002

# Run 1: SRR21832358
echo '  Downloading SRR21832358...'
prefetch SRR21832358 || { echo 'Failed to prefetch SRR21832358'; exit 1; }
fasterq-dump SRR21832358 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832358'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832358_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832358_3.fastq primary_tumor_002_S1_L001_R1_001.fastq
    mv SRR21832358_4.fastq primary_tumor_002_S1_L001_R2_001.fastq
    rm -f SRR21832358_1.fastq SRR21832358_2.fastq
elif [ -f SRR21832358_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832358_2.fastq primary_tumor_002_S1_L001_R1_001.fastq
    mv SRR21832358_3.fastq primary_tumor_002_S1_L001_R2_001.fastq
    rm -f SRR21832358_1.fastq
elif [ -f SRR21832358_1.fastq ] && [ -f SRR21832358_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832358_1.fastq primary_tumor_002_S1_L001_R1_001.fastq
    mv SRR21832358_2.fastq primary_tumor_002_S1_L001_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832358'
    ls -lh SRR21832358*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_002_S1_L001_R1_001.fastq
gzip primary_tumor_002_S1_L001_R2_001.fastq

rm -rf SRR21832358

# Run 2: SRR21832359
echo '  Downloading SRR21832359...'
prefetch SRR21832359 || { echo 'Failed to prefetch SRR21832359'; exit 1; }
fasterq-dump SRR21832359 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832359'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832359_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832359_3.fastq primary_tumor_002_S1_L002_R1_001.fastq
    mv SRR21832359_4.fastq primary_tumor_002_S1_L002_R2_001.fastq
    rm -f SRR21832359_1.fastq SRR21832359_2.fastq
elif [ -f SRR21832359_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832359_2.fastq primary_tumor_002_S1_L002_R1_001.fastq
    mv SRR21832359_3.fastq primary_tumor_002_S1_L002_R2_001.fastq
    rm -f SRR21832359_1.fastq
elif [ -f SRR21832359_1.fastq ] && [ -f SRR21832359_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832359_1.fastq primary_tumor_002_S1_L002_R1_001.fastq
    mv SRR21832359_2.fastq primary_tumor_002_S1_L002_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832359'
    ls -lh SRR21832359*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_002_S1_L002_R1_001.fastq
gzip primary_tumor_002_S1_L002_R2_001.fastq

rm -rf SRR21832359

cd ../../..
echo 'Completed primary_tumor_002'
echo ''


# ============================================================
# Sample: primary_tumor_003 (Original: GSM6619237)
# Tumor type: primary tumor, Age: 67 years
# ============================================================

echo 'Processing primary_tumor_003...'
mkdir -p fastq/cellranger/primary_tumor_003
cd fastq/cellranger/primary_tumor_003

# Run 1: SRR21832356
echo '  Downloading SRR21832356...'
prefetch SRR21832356 || { echo 'Failed to prefetch SRR21832356'; exit 1; }
fasterq-dump SRR21832356 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832356'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832356_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832356_3.fastq primary_tumor_003_S1_L001_R1_001.fastq
    mv SRR21832356_4.fastq primary_tumor_003_S1_L001_R2_001.fastq
    rm -f SRR21832356_1.fastq SRR21832356_2.fastq
elif [ -f SRR21832356_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832356_2.fastq primary_tumor_003_S1_L001_R1_001.fastq
    mv SRR21832356_3.fastq primary_tumor_003_S1_L001_R2_001.fastq
    rm -f SRR21832356_1.fastq
elif [ -f SRR21832356_1.fastq ] && [ -f SRR21832356_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832356_1.fastq primary_tumor_003_S1_L001_R1_001.fastq
    mv SRR21832356_2.fastq primary_tumor_003_S1_L001_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832356'
    ls -lh SRR21832356*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_003_S1_L001_R1_001.fastq
gzip primary_tumor_003_S1_L001_R2_001.fastq

rm -rf SRR21832356

# Run 2: SRR21832357
echo '  Downloading SRR21832357...'
prefetch SRR21832357 || { echo 'Failed to prefetch SRR21832357'; exit 1; }
fasterq-dump SRR21832357 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832357'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832357_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832357_3.fastq primary_tumor_003_S1_L002_R1_001.fastq
    mv SRR21832357_4.fastq primary_tumor_003_S1_L002_R2_001.fastq
    rm -f SRR21832357_1.fastq SRR21832357_2.fastq
elif [ -f SRR21832357_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832357_2.fastq primary_tumor_003_S1_L002_R1_001.fastq
    mv SRR21832357_3.fastq primary_tumor_003_S1_L002_R2_001.fastq
    rm -f SRR21832357_1.fastq
elif [ -f SRR21832357_1.fastq ] && [ -f SRR21832357_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832357_1.fastq primary_tumor_003_S1_L002_R1_001.fastq
    mv SRR21832357_2.fastq primary_tumor_003_S1_L002_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832357'
    ls -lh SRR21832357*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_003_S1_L002_R1_001.fastq
gzip primary_tumor_003_S1_L002_R2_001.fastq

rm -rf SRR21832357

# Run 3: SRR21832364
echo '  Downloading SRR21832364...'
prefetch SRR21832364 || { echo 'Failed to prefetch SRR21832364'; exit 1; }
fasterq-dump SRR21832364 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832364'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832364_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832364_3.fastq primary_tumor_003_S1_L003_R1_001.fastq
    mv SRR21832364_4.fastq primary_tumor_003_S1_L003_R2_001.fastq
    rm -f SRR21832364_1.fastq SRR21832364_2.fastq
elif [ -f SRR21832364_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832364_2.fastq primary_tumor_003_S1_L003_R1_001.fastq
    mv SRR21832364_3.fastq primary_tumor_003_S1_L003_R2_001.fastq
    rm -f SRR21832364_1.fastq
elif [ -f SRR21832364_1.fastq ] && [ -f SRR21832364_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832364_1.fastq primary_tumor_003_S1_L003_R1_001.fastq
    mv SRR21832364_2.fastq primary_tumor_003_S1_L003_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832364'
    ls -lh SRR21832364*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_003_S1_L003_R1_001.fastq
gzip primary_tumor_003_S1_L003_R2_001.fastq

rm -rf SRR21832364

# Run 4: SRR21832369
echo '  Downloading SRR21832369...'
prefetch SRR21832369 || { echo 'Failed to prefetch SRR21832369'; exit 1; }
fasterq-dump SRR21832369 --split-files --include-technical --threads 8 || { echo 'Failed to dump SRR21832369'; exit 1; }

# Detect read layout and rename to Cell Ranger format:
# 4 files (_1,_2,_3,_4): dual-index 10x  -> _3=R1 (barcode+UMI), _4=R2 (cDNA)
# 3 files (_1,_2,_3):    single-index 10x -> _2=R1 (barcode+UMI), _3=R2 (cDNA)
# 2 files (_1,_2):        standard paired  -> _1=R1, _2=R2
if [ -f SRR21832369_4.fastq ]; then
    echo '  4-file layout: using _3=R1, _4=R2'
    mv SRR21832369_3.fastq primary_tumor_003_S1_L004_R1_001.fastq
    mv SRR21832369_4.fastq primary_tumor_003_S1_L004_R2_001.fastq
    rm -f SRR21832369_1.fastq SRR21832369_2.fastq
elif [ -f SRR21832369_3.fastq ]; then
    echo '  3-file layout: using _2=R1, _3=R2'
    mv SRR21832369_2.fastq primary_tumor_003_S1_L004_R1_001.fastq
    mv SRR21832369_3.fastq primary_tumor_003_S1_L004_R2_001.fastq
    rm -f SRR21832369_1.fastq
elif [ -f SRR21832369_1.fastq ] && [ -f SRR21832369_2.fastq ]; then
    echo '  2-file layout: using _1=R1, _2=R2'
    mv SRR21832369_1.fastq primary_tumor_003_S1_L004_R1_001.fastq
    mv SRR21832369_2.fastq primary_tumor_003_S1_L004_R2_001.fastq
else
    echo 'ERROR: Cannot find expected FASTQ files for SRR21832369'
    ls -lh SRR21832369*.fastq 2>/dev/null || echo 'No fastq files found'
    exit 1
fi

echo '  Compressing...'
gzip primary_tumor_003_S1_L004_R1_001.fastq
gzip primary_tumor_003_S1_L004_R2_001.fastq

rm -rf SRR21832369

cd ../../..
echo 'Completed primary_tumor_003'
echo ''


echo 'All samples downloaded and renamed!'
echo 'Directory structure:'
tree fastq/ -L 2
