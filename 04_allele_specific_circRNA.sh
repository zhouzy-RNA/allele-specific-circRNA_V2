### Preparatory work ###
#1. Exons in the GTF file
#2. The autosomal genome is converted into one line
#3. N-masked autosomal genome
#4. SNP information of each autosomal
		#1	169	170	rs1112347109	G	T
		#1	369	370	rs1108162481	G	A
		#1	481	482	rs1107705765	T	C
		#1	528	529	.	G	A
		#1	554	555	rs1113680398	C	G
		#1	1081	1082	rs1111334233	G	T
#5. Autosomal number





##### What this script contains? #####
#1. The sequences of the reads aligned to RNAs were segmented according to genomic location.
#2. The sequences of the reads aligned to RNAs were converted to genomic location to overlop with SNPs
#3. The segmented read intersect with SNP to distinguish different alleles
#4. The number of reads of different alleles was at least 3
#5. ### Allelic specificity was validated by odds ratio test





sample="DKM2DNF5_4"
chr_number=18
chr_seq="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sscrofa11.1.dna.chr.fa"
gtf="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sus_scrofa.Sscrofa11.1.90.chromsome.gtf"
N_masked_ref_seq="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sscrofa11.1.dna.chr.masked.fa"

vcf_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/vcf"
script_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/script"
data_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/data"

bwa_path="/opt/service/fermi.kit"
samtools_path="/opt/service/miniconda3/envs/rex_env/bin"
Rscript_pasth="/opt/service/miniconda3/envs/rex_env/bin"
bedtools_path="/opt/service/miniconda3/envs/rex_env/bin"
CIRI2_path="/data_group/zhouzhongyin/zhouzhongyin1/software/CIRI_v2.0.6"

out_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/out"
ciri2_out_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ciri2_out"

	cd ${out_path}
	cd $sample



	### The sequences of the read1 aligned to circRNAs were segmented according to genomic location (reverse strand)
	cut -f2 ${sample}.circ2.info.intersect_exon.neg.bed.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]{if ($(NF-1)=="junction"){printf $1"\t"$2"\t";for (i=3;i<=NF-4;i+=2){printf $i"\t"$(i+1)"\t"};print $(NF-3)"\t"$(NF-2)}else{print $0}}' - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg ${sample}.1.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads

	cut -f3 ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads | rev | sed "s/A/#/g;s/G/%/g;s/C/G/g;s/T/A/g;s/#/T/g;s/%/C/g" > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.temp_rev

	cut -f1,2 ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads | paste - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.temp_rev > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=1;for (i=3;i<=NF;i+=2){print $1"\t"$2"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq





	### The sequences of the read1 aligned to circRNAs were segmented according to genomic location (forword strand)
	cut -f2 ${sample}.circ2.info.intersect_exon.pos.bed.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]{if ($(NF-1)=="junction"){printf $1"\t"$2"\t";for (i=3;i<=NF-4;i+=2){printf $i"\t"$(i+1)"\t"};print $(NF-3)"\t"$(NF-2)}else{print $0}}' - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos ${sample}.1.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=1;for (i=3;i<=NF;i+=2){print $1"\t"$2"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq





	### The sequences of the read2 aligned to circRNAs were segmented according to genomic location (reverse strand)
	cut -f2 ${sample}.circ2.info.intersect_exon.neg.bed.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]{if ($(NF-1)=="junction"){printf $1"\t"$2"\t";for (i=3;i<=NF-4;i+=2){printf $i"\t"$(i+1)"\t"};print $(NF-3)"\t"$(NF-2)}else{print $0}}'  - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg ${sample}.2.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads

	cut -f3 ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads | rev | sed "s/A/#/g;s/G/%/g;s/C/G/g;s/T/A/g;s/#/T/g;s/%/C/g" > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.temp_rev

	cut -f1,2 ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads | paste - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.temp_rev > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=1;for (i=3;i<=NF;i+=2){print $1"\t"$2"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq





	### The sequences of the read2 aligned to circRNAs were segmented according to genomic location (forword strand)
	cut -f2 ${sample}.circ2.info.intersect_exon.pos.bed.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]{if ($(NF-1)=="junction"){printf $1"\t"$2"\t";for (i=3;i<=NF-4;i+=2){printf $i"\t"$(i+1)"\t"};print $(NF-3)"\t"$(NF-2)}else{print $0}}' - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos ${sample}.2.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=1;for (i=3;i<=NF;i+=2){print $1"\t"$2"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq





	### The sequences of the read1 aligned to linear RNA were segmented according to genomic location (reverse strand)
	cut -f2 ${sample}.circ2.linear_RNA.bed.neg.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]&&a[$2]&&$3!="left"&&$3!="right"{print $0}' - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg ${sample}.1.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads

	cut -f3 ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads | rev | sed "s/A/#/g;s/G/%/g;s/C/G/g;s/T/A/g;s/#/T/g;s/%/C/g" > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev

	cut -f1,2 ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads | paste - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=$4;for (i=5;i<=NF-2;i+=2){print $1"\t"$NF"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq





	### The sequences of the read1 aligned to linear RNA were segmented according to genomic location (forword strand)
	cut -f2 ${sample}.circ2.linear_RNA.bed.pos.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]&&$3!="left"&&$3!="right"{print $0}' - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos ${sample}.1.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=$4;for (i=5;i<=NF-2;i+=2){print $1"\t"$NF"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq





	### The sequences of the read2 aligned to linear RNA were segmented according to genomic location (reverse strand)
	cut -f2 ${sample}.circ2.linear_RNA.bed.neg.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]&&$3!="left"&&$3!="right"{print $0}' - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg ${sample}.2.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads

	cut -f3 ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads | rev | sed "s/A/#/g;s/G/%/g;s/C/G/g;s/T/A/g;s/#/T/g;s/%/C/g" > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev

	cut -f1,2 ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads | paste - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=$4;for (i=5;i<=NF-2;i+=2){print $1"\t"$NF"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq





	### The sequences of the read2 aligned to linear RNA were segmented according to genomic location (forword strand)
	cut -f2 ${sample}.circ2.linear_RNA.bed.pos.one_line | sort |uniq | sed "s/_[0-9]*//g" | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$2]&&$3!="left"&&$3!="right"{print $0}' - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos

	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1"\t"$3"\t"$6"\t"$10}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos ${sample}.2.all.sorted.sam.temp_noheader | sed "s/_[0-9]*//g" | awk '{a[$1]++;if (a[$1]==1){c[$1]=$4;print $0}else{print $1"\t"$2"\t"$3"\t"c[$1]}}' | awk '$3!~/[*IDSHP=X]/{print}' | cut -f1,2,4 > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads

	awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{p=$4;for (i=5;i<=NF-2;i+=2){print $1"\t"$NF"\t"$i"\t"$(i+1)"\t"substr(a[$1"\t"$2],p,$(i+1)-$i+1);p+=$(i+1)-$i+1}}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq

	rm ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.reads.rev ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.reads ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos  ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.temp_rev ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.reads.rev   ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos  ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.reads


	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq

	cat ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq

	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq

	cat ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq

	rm ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.neg.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.pos.seq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.neg.seq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.pos.seq





	### The sequences of the reads aligned to RNAs were converted to genomic location to overlop with SNPs
	awk '{match($2,"([0-9]+):",arr);print arr[1]"\t"$3"\t"$4"\t"$2"\t"$1"\t"$5}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq.bed

	awk '{match($2,"([0-9]+):",arr);print arr[1]"\t"$3"\t"$4"\t"$2"\t"$1"\t"$5}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq.bed

	awk '{match($2,"([0-9]+):",arr);print arr[1]"\t"$3"\t"$4"\t"$2"\t"$1"\t"$5}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq.bed

	awk '{match($2,"([0-9]+):",arr);print arr[1]"\t"$3"\t"$4"\t"$2"\t"$1"\t"$5}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq.bed






	### The segmented read intersect with SNP to distinguish different alleles
	rm ${sample}.1.circular_RNA.intersect.vcf.all_chr.bed ${sample}.2.circular_RNA.intersect.vcf.all_chr.bed ${sample}.1.linear_RNA.intersect.vcf.all_chr.bed ${sample}.2.linear_RNA.intersect.vcf.all_chr.bed
	for ((i=1;i<=$chr_number;i++))
	do

		awk '$1==a{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}' a=$i ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq.bed | ${bedtools_path}/bedtools intersect -wa -wb -a - -b ${vcf_path}/chr${i}_overlapping.chr_posit.vcf.bed > ${sample}.1.circular_RNA.intersect.vcf.bed
		cat ${sample}.1.circular_RNA.intersect.vcf.bed >> ${sample}.1.circular_RNA.intersect.vcf.all_chr.bed
		rm ${sample}.1.circular_RNA.intersect.vcf.bed

		awk '$1==a{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}' a=$i ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.seq.bed | ${bedtools_path}/bedtools intersect -wa -wb -a - -b ${vcf_path}/chr${i}_overlapping.chr_posit.vcf.bed > ${sample}.2.circular_RNA.intersect.vcf.bed
		cat ${sample}.2.circular_RNA.intersect.vcf.bed >> ${sample}.2.circular_RNA.intersect.vcf.all_chr.bed
		rm ${sample}.2.circular_RNA.intersect.vcf.bed

		awk '$1==a{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}' a=$i ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq.bed | ${bedtools_path}/bedtools intersect -wa -wb -a - -b ${vcf_path}/chr${i}_overlapping.chr_posit.vcf.bed > ${sample}.1.linear_RNA.intersect.vcf.bed
		cat ${sample}.1.linear_RNA.intersect.vcf.bed >> ${sample}.1.linear_RNA.intersect.vcf.all_chr.bed
		rm ${sample}.1.linear_RNA.intersect.vcf.bed

		awk '$1==a{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}' a=$i ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region.seq.bed | ${bedtools_path}/bedtools intersect -wa -wb -a - -b ${vcf_path}/chr${i}_overlapping.chr_posit.vcf.bed > ${sample}.2.linear_RNA.intersect.vcf.bed
		cat ${sample}.2.linear_RNA.intersect.vcf.bed >> ${sample}.2.linear_RNA.intersect.vcf.all_chr.bed
		rm ${sample}.2.linear_RNA.intersect.vcf.bed

	done

	cat ${sample}.1.circular_RNA.intersect.vcf.all_chr.bed ${sample}.2.circular_RNA.intersect.vcf.all_chr.bed | awk '{$2=$2+1;print $0}' | awk '{print $4"\t"$5"\t"substr($6,$9-$2+1,1)"\t"$9"\t"$10"\t"$11"\t"$12}' | sort |uniq | awk '$3==$6||$3==$7{print}' | cut -f1,3- | awk '{a[$0]++}END{for (i in a)print i"\t"a[i]}' | sort -k1,1 -k3,3n > ${sample}.circular_RNA.out

	cat ${sample}.1.linear_RNA.intersect.vcf.all_chr.bed ${sample}.2.linear_RNA.intersect.vcf.all_chr.bed | awk '{$2=$2+1;print $0}' | awk '{print $4"\t"$5"\t"substr($6,$9-$2+1,1)"\t"$9"\t"$10"\t"$11"\t"$12}' | sort |uniq | awk '$3==$6||$3==$7{print}' | cut -f1,3- | awk '{a[$0]++}END{for (i in a)print i"\t"a[i]}' | sort -k1,1 -k3,3n > ${sample}.linear_RNA.out





	### The number of reads of different alleles was at least 3
	awk 'NR==FNR&&$NF>=3{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]==2{print}' ${sample}.circular_RNA.out ${sample}.circular_RNA.out > ${sample}.circular_RNA.out.3.need

	awk 'NR==FNR&&$NF>=3{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]==2{print}' ${sample}.linear_RNA.out ${sample}.linear_RNA.out > ${sample}.linear_RNA.out.3.need

	awk 'NR==FNR{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]{print}' ${sample}.circular_RNA.out.3.need ${sample}.linear_RNA.out.3.need | awk 'NR==FNR{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]{print $0}' - ${sample}.circular_RNA.out.3.need > ${sample}.circular_RNA.out.3.need.with_linear

	awk 'NR==FNR{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]{print}' ${sample}.circular_RNA.out.3.need ${sample}.linear_RNA.out.3.need | awk 'NR==FNR{a[$1"\t"$3]++}NR>FNR&&a[$1"\t"$3]{print $0}' - ${sample}.linear_RNA.out.3.need > ${sample}.linear_RNA.out.3.need.with_linear

	awk '{if (FNR%2==1){temp=$7}else{if ($2==$5){print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"temp}else{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"temp"\t"$7}}}' ${sample}.circular_RNA.out.3.need.with_linear > ${sample}.circular_RNA.out.3.need.with_linear.ref_alt

	awk '{if (FNR%2==1){temp=$7}else{if ($2==$5){print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"temp}else{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"temp"\t"$7}}}' ${sample}.linear_RNA.out.3.need.with_linear > ${sample}.linear_RNA.out.3.need.with_linear.ref_alt

	awk 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4"\t"$5]=$6"\t"$7}NR>FNR{print $0"\t"a[$1"\t"$2"\t"$3"\t"$4"\t"$5]}' ${sample}.circular_RNA.out.3.need.with_linear.ref_alt ${sample}.linear_RNA.out.3.need.with_linear.ref_alt > ${sample}.linear_circular_RNA.out.3.need.ref_alt

	rm ${sample}.1.all.sorted.sam.temp_noheader ${sample}.2.all.sorted.sam.temp_noheader ${sample}.1.all.sorted.sam.temp_noheader.complete_map ${sample}.2.all.sorted.sam.temp_noheader.complete_map ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need





	### Allelic specificity was validated by odds ratio test
	${Rscript_pasth}/Rscript ${script_path}/oddsratio.R ${sample}


	grep "|" ${sample}.linear_circular_RNA.out.3.need.ref_alt.oddsratio | awk '$10>1&&$11>1&&$12>1{print}' | sed "s/\//\t/g;s/oddsratio:/\t/g" | awk '{print a"\t"$0}' a=${sample} - > ${sample}.linear_circular_RNA.out.3.need.ref_alt.oddsratio.ASE-circRNA
