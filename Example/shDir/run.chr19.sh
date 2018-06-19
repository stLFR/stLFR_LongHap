perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:1-10000000 300000 1 1 outDir_tmp/chr19_1_10000000.log 2>outDir_tmp/chr19_1_10000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:9500001-20000000 300000 1 1 outDir_tmp/chr19_9500001_20000000.log 2>outDir_tmp/chr19_9500001_20000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:19500001-30000000 300000 1 1 outDir_tmp/chr19_19500001_30000000.log 2>outDir_tmp/chr19_19500001_30000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:29500001-40000000 300000 1 1 outDir_tmp/chr19_29500001_40000000.log 2>outDir_tmp/chr19_29500001_40000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:39500001-50000000 300000 1 1 outDir_tmp/chr19_39500001_50000000.log 2>outDir_tmp/chr19_39500001_50000000.error
perl /home-02/ysun/LongHap/Bin/LongHap.pl input.bam input.vcf chr19:49500001-59128983 300000 1 1 outDir_tmp/chr19_49500001_59128983.log 2>outDir_tmp/chr19_49500001_59128983.error
perl /home-02/ysun/LongHap/Bin/merge.pl outDir_tmp chr19 10000000 outDir_final/chr19.log
