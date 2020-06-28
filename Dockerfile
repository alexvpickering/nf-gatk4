FROM broadinstitute/gatk:4.1.7.0

# hg38 with alt contigs o Homo_sapiens_assembly38.fasta.64.alt is downloaded 
RUN mkdir ~/hg38 && cd ~/hg38 \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt \
    && wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

RUN apt-get update && apt-get install bwa=0.7.12-5

RUN cd ~/hg38 \
    && wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

RUN cd ~/hg38 \
    && wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
    && wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi \
    && wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/hapmap_3.3.hg38.vcf.gz \
    && wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
