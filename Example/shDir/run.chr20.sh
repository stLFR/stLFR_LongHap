perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:1-10000000 300000 1 1 outDir_tmp/chr20_1_10000000.log 2>outDir_tmp/chr20_1_10000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:9500001-20000000 300000 1 1 outDir_tmp/chr20_9500001_20000000.log 2>outDir_tmp/chr20_9500001_20000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:19500001-30000000 300000 1 1 outDir_tmp/chr20_19500001_30000000.log 2>outDir_tmp/chr20_19500001_30000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:29500001-40000000 300000 1 1 outDir_tmp/chr20_29500001_40000000.log 2>outDir_tmp/chr20_29500001_40000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:39500001-50000000 300000 1 1 outDir_tmp/chr20_39500001_50000000.log 2>outDir_tmp/chr20_39500001_50000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:49500001-60000000 300000 1 1 outDir_tmp/chr20_49500001_60000000.log 2>outDir_tmp/chr20_49500001_60000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr20:59500001-63025520 300000 1 1 outDir_tmp/chr20_59500001_63025520.log 2>outDir_tmp/chr20_59500001_63025520.error
perl /home-02/ysun/LongHap/Bin/merge.pl outDir_tmp chr20 10000000 outDir_final/chr20.log
