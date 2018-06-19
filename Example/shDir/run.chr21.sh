perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr21:1-10000000 300000 1 1 outDir_tmp/chr21_1_10000000.log 2>outDir_tmp/chr21_1_10000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr21:9500001-20000000 300000 1 1 outDir_tmp/chr21_9500001_20000000.log 2>outDir_tmp/chr21_9500001_20000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr21:19500001-30000000 300000 1 1 outDir_tmp/chr21_19500001_30000000.log 2>outDir_tmp/chr21_19500001_30000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr21:29500001-40000000 300000 1 1 outDir_tmp/chr21_29500001_40000000.log 2>outDir_tmp/chr21_29500001_40000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr21:39500001-48129895 300000 1 1 outDir_tmp/chr21_39500001_48129895.log 2>outDir_tmp/chr21_39500001_48129895.error
perl /home-02/ysun/LongHap/Bin/merge.pl outDir_tmp chr21 10000000 outDir_final/chr21.log
