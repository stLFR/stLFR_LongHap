perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:1-10000000 300000 1 1 outDir_tmp/chr22_1_10000000.log 2>outDir_tmp/chr22_1_10000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:9500001-20000000 300000 1 1 outDir_tmp/chr22_9500001_20000000.log 2>outDir_tmp/chr22_9500001_20000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:19500001-30000000 300000 1 1 outDir_tmp/chr22_19500001_30000000.log 2>outDir_tmp/chr22_19500001_30000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:29500001-40000000 300000 1 1 outDir_tmp/chr22_29500001_40000000.log 2>outDir_tmp/chr22_29500001_40000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:39500001-50000000 300000 1 1 outDir_tmp/chr22_39500001_50000000.log 2>outDir_tmp/chr22_39500001_50000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr22:49500001-51304566 300000 1 1 outDir_tmp/chr22_49500001_51304566.log 2>outDir_tmp/chr22_49500001_51304566.error
perl /home-02/ysun/LongHap/Bin/merge.pl outDir_tmp chr22 10000000 outDir_final/chr22.log
