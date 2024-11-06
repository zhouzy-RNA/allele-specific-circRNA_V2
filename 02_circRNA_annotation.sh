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
#1. circRNA annotation by CIRI2
#2. Exhaustive transcripts, exons of which cover both back-splcing sites. Exons of these transcripts between back-splcing sites were candidate of structure of circular RNA 
#3. The uniq structure of circular RNA were preserved. Exons of different transcripts between back-splcing sites may be the same
#4. the sequence of circular RNA and linear RNA were Extracted from N-masked genome to generated reference transcriptome
#5. Build the reference transcriptome index





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
	mkdir $sample
	cd $sample

	### circRNA annotation
	${bwa_path}/bwa mem -T 19 -t 5 \
		$chr_seq \
		${data_path}/${sample}_1.fq.gz \
		${data_path}/${sample}_2.fq.gz \
		> ${sample}.sam

	perl ${CIRI2_path}/CIRI2.pl \
		-I ${sample}.sam -O ${ciri2_out_path}/${sample}.ciri2.out \
		-F $chr_seq \
		-A $gtf

	rm ${sample}.sam

	### Genomic location and gene information of exon circular RNA
	awk '$9=="exon"{print}' ${ciri2_out_path}/${sample}.ciri2.out | awk '{print $2"\t"$3"\t"$4"\t"$11"\t"$1"\t"$10}' | sed "s/,//g" | grep "^[1-9]" > ${sample}.circ2.info
	#7       26573830        26664601        -       7:26573830|26664601     ENSSSCG00000001485





	### Exon information in circular RNA (reverse strand)
	${bedtools_path}/bedtools intersect -wa -wb -a ${sample}.circ2.info -b ${gtf}.exon.bed | awk '$6==$11{print $0}' | awk '$4=="-"{print}' | sort -k5,5 -k12,12 -k13,13nr > ${sample}.circ2.info.intersect_exon.neg
	#7       26573830        26664601        -       7:26573830|26664601     ENSSSCG00000001485      7       26599359        26599601        -       ENSSSCG00000001485      ENSSSCT00000001657      11

	### Exhaustive transcripts, these exons of transcripts cover both back-splcing sites
	### Exon information of these transcripts between back-splcing sites were preserved
	awk 'BEGIN{m=0;n=0}NR==FNR{a[$5"\t"$12]++}NR>FNR{b[$5"\t"$12]++;
		if (b[$5"\t"$12]==1){
			if ($2>=$8&&$2<=$9){
				m=1
			}
		};
		if (b[$5"\t"$12]==a[$5"\t"$12]){
			if ($3>=$8&&$3<=$9){
				n=1
			}
			if (m==1&&n==1){
				print $5"\t"$12
			}
			m=0;n=0
		}
	}' ${sample}.circ2.info.intersect_exon.neg ${sample}.circ2.info.intersect_exon.neg | awk 'NR==FNR{a[$0]++}NR>FNR&&a[$5"\t"$12]{print $0}' - ${sample}.circ2.info.intersect_exon.neg | awk 'NR==FNR{a[$5"\t"$12]++}NR>FNR{b[$5"\t"$12]++;if (b[$5"\t"$12]==1){print $7"\t"$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}else if (b[$5"\t"$12]==a[$5"\t"$12]){print $7"\t"$8"\t"$3"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}else{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}}' ${sample}.circ2.info.intersect_exon.neg - | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$7"\t"$6"\t"$5}' > ${sample}.circ2.info.intersect_exon.neg.temp_bed
	#7       26573830        26573883        -       7:26573830|26664601     10      ENSSSCT00000029860      ENSSSCG00000001485

	### The circular RNA exons were sorted and numbered according to the ensemble annotation
	sort -k5,5 -k7,7 -k6,6n ${sample}.circ2.info.intersect_exon.neg.temp_bed | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){b[$5"\t"$7]=$6};$6=$6-b[$5"\t"$7]+1;print $0}' | sed "s/ /\t/g" > ${sample}.circ2.info.intersect_exon.neg.temp_bed2

	### The exon information for each transcript is converted into a single line
	awk '{if (a[$5"\t"$7]){a[$5"\t"$7]=a[$5"\t"$7]"\t"$2"\t"$3}else{a[$5"\t"$7]=$2"\t"$3}}END{for (i in a){print i"\t"a[i]}}' ${sample}.circ2.info.intersect_exon.neg.temp_bed2 | sort > ${sample}.circ2.info.intersect_exon.neg.temp_bed3

	### The uniq gene structures were preserved
	cut -f3- ${sample}.circ2.info.intersect_exon.neg.temp_bed3 | awk '{a[$0]++;print a[$0]"\t"$0}' | awk 'NR==FNR{if ($1==1){a[FNR]++}}NR>FNR&&a[FNR]{print}' - ${sample}.circ2.info.intersect_exon.neg.temp_bed3 > ${sample}.circ2.info.intersect_exon.neg.temp_bed4

	### The annotation of circular RNA was generated from the uniq gene structure
	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$5"\t"$7]{print $0}' ${sample}.circ2.info.intersect_exon.neg.temp_bed4 ${sample}.circ2.info.intersect_exon.neg.temp_bed2 | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){b[$5]++};$5=$5"_"b[$5];print $0}' | sed "s/ /\t/g" > ${sample}.circ2.info.intersect_exon.neg.bed

	rm ${sample}.circ2.info.intersect_exon.neg.temp_bed ${sample}.circ2.info.intersect_exon.neg.temp_bed2 ${sample}.circ2.info.intersect_exon.neg.temp_bed3 ${sample}.circ2.info.intersect_exon.neg.temp_bed4





	### Exon information in circular RNA (forword strand)
	${bedtools_path}/bedtools intersect -wa -wb -a ${sample}.circ2.info -b ${gtf}.exon.bed | awk '$6==$11{print $0}' | awk '$4=="+"{print}' | sort -k5,5 -k12,12 -k13,13n > ${sample}.circ2.info.intersect_exon.pos

	### Exhaustive transcripts, these exons of transcripts cover both back-splcing sites
	### Exon information of these transcripts between back-splcing sites were preserved
	awk 'BEGIN{m=0;n=0}NR==FNR{a[$5"\t"$12]++}NR>FNR{b[$5"\t"$12]++;
		if (b[$5"\t"$12]==1){
			if ($2>=$8&&$2<=$9){
				m=1
			}
		};
		if (b[$5"\t"$12]==a[$5"\t"$12]){
			if ($3>=$8&&$3<=$9){
				n=1
			}
			if (m==1&&n==1){
				print $5"\t"$12
			}
			m=0;n=0
		}
	}' ${sample}.circ2.info.intersect_exon.pos ${sample}.circ2.info.intersect_exon.pos | awk 'NR==FNR{a[$0]++}NR>FNR&&a[$5"\t"$12]{print $0}' - ${sample}.circ2.info.intersect_exon.pos | awk 'NR==FNR{a[$5"\t"$12]++}NR>FNR{b[$5"\t"$12]++;if (b[$5"\t"$12]==1){print $7"\t"$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}else if (b[$5"\t"$12]==a[$5"\t"$12]){print $7"\t"$8"\t"$3"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}else{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$5}}' ${sample}.circ2.info.intersect_exon.pos - | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$7"\t"$6"\t"$5}' > ${sample}.circ2.info.intersect_exon.pos.temp_bed

	### The circular RNA exons were sorted and numbered according to the ensemble annotation
	sort -k5,5 -k7,7 -k6,6n ${sample}.circ2.info.intersect_exon.pos.temp_bed | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){b[$5"\t"$7]=$6};$6=$6-b[$5"\t"$7]+1;print $0}' | sed "s/ /\t/g" > ${sample}.circ2.info.intersect_exon.pos.temp_bed2

	### The exon information for each transcript is converted into a single line
	awk '{if (a[$5"\t"$7]){a[$5"\t"$7]=a[$5"\t"$7]"\t"$2"\t"$3}else{a[$5"\t"$7]=$2"\t"$3}}END{for (i in a){print i"\t"a[i]}}' ${sample}.circ2.info.intersect_exon.pos.temp_bed2 | sort > ${sample}.circ2.info.intersect_exon.pos.temp_bed3

	### The uniq gene structures were preserved
	cut -f3- ${sample}.circ2.info.intersect_exon.pos.temp_bed3 | awk '{a[$0]++;print a[$0]"\t"$0}' | awk 'NR==FNR{if ($1==1){a[FNR]++}}NR>FNR&&a[FNR]{print}' - ${sample}.circ2.info.intersect_exon.pos.temp_bed3 > ${sample}.circ2.info.intersect_exon.pos.temp_bed4

	### The annotation of circular RNA was generated from the uniq gene structure
	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$5"\t"$7]{print $0}' ${sample}.circ2.info.intersect_exon.pos.temp_bed4 ${sample}.circ2.info.intersect_exon.pos.temp_bed2 | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){b[$5]++};$5=$5"_"b[$5];print $0}' | sed "s/ /\t/g" > ${sample}.circ2.info.intersect_exon.pos.bed

	rm ${sample}.circ2.info.intersect_exon.pos.temp_bed ${sample}.circ2.info.intersect_exon.pos.temp_bed2 ${sample}.circ2.info.intersect_exon.pos.temp_bed3 ${sample}.circ2.info.intersect_exon.pos.temp_bed4






	### Combining annotated information of reverse and forword strand circRNA
	cat ${sample}.circ2.info.intersect_exon.neg.bed ${sample}.circ2.info.intersect_exon.pos.bed > ${sample}.circ2.info.intersect_exon.bed






	### Linear transcript structure
	awk 'NR==FNR{a[$6]++}NR>FNR&&a[$5]{print $0}' ${sample}.circ2.info ${gtf}.exon.bed > ${sample}.circ2.linear_RNA.temp_bed

	awk '{a[$5"\t"$6]++;if (a[$5"\t"$6]==1){b[$5]++};$5=$5"_"b[$5];print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$6}' ${sample}.circ2.linear_RNA.temp_bed > ${sample}.circ2.linear_RNA.bed






	### Extract circRNA sequence
	perl ${script_path}/ExtractTranscriptSequence.pl \
		-genome ${N_masked_ref_seq} \
		-transcript ${sample}.circ2.info.intersect_exon.bed \
		-output ${sample}.circ2.info.intersect_exon.bed.temp

	perl ${script_path}/ConcatExons.pl \
		-input_fasta ${sample}.circ2.info.intersect_exon.bed.temp \
		-output  ${sample}.circ2.info.intersect_exon.bed.temp2
	
	### repeat circRNA sequence
	awk '{if (FNR%2==1){print $1}else{print $1""$1}}' ${sample}.circ2.info.intersect_exon.bed.temp2 > ${sample}.circ2.info.intersect_exon.bed.fa

	rm ${sample}.circ2.info.intersect_exon.bed.temp ${sample}.circ2.info.intersect_exon.bed.temp2






	### Extract linear RNA sequence
	perl ${script_path}/ExtractTranscriptSequence.pl \
		-genome ${N_masked_ref_seq} \
		-transcript ${sample}.circ2.linear_RNA.bed \
		-output ${sample}.circ2.linear_RNA.bed.temp

	perl ${script_path}/ConcatExons.pl \
		-input_fasta ${sample}.circ2.linear_RNA.bed.temp \
		-output ${sample}.circ2.linear_RNA.bed.fa

	rm ${sample}.circ2.linear_RNA.bed.temp






	### Combining the sequence of circular RNA and linear RNA
	cat ${sample}.circ2.info.intersect_exon.bed.fa ${sample}.circ2.linear_RNA.bed.fa > ${sample}.circRNA_linear_RNA.fa

	### Build the reference transcriptome index
	${bwa_path}/bwa index -a is ${sample}.circRNA_linear_RNA.fa






	### The exon information for circular RNA and linear RNA is arranged in a line
	awk '$4=="-"{print}' ${sample}.circ2.linear_RNA.bed | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){if (FNR==1){printf $7"\t"$5"\t"$1"\t"$4"\t"$3"\t"$2}else{print "";printf $7"\t"$5"\t"$1"\t"$4"\t"$3"\t"$2}}else{printf "\t"$3"\t"$2}}END{print ""}' > ${sample}.circ2.linear_RNA.bed.neg.one_line

	awk '$4=="+"{print}' ${sample}.circ2.linear_RNA.bed | awk '{a[$5"\t"$7]++;if (a[$5"\t"$7]==1){if (FNR==1){printf $7"\t"$5"\t"$1"\t"$4"\t"$2"\t"$3}else{print "";printf $7"\t"$5"\t"$1"\t"$4"\t"$2"\t"$3}}else{printf "\t"$2"\t"$3}}END{print ""}' > ${sample}.circ2.linear_RNA.bed.pos.one_line

	awk '{a[$7"\t"$5]++;if (a[$7"\t"$5]==1){if (FNR==1){printf $7"\t"$5"\t"$1"\t"$4"\t"$3"\t"$2}else{print "";printf $7"\t"$5"\t"$1"\t"$4"\t"$3"\t"$2}}else{printf "\t"$3"\t"$2}}END{print ""}' ${sample}.circ2.info.intersect_exon.neg.bed > ${sample}.circ2.info.intersect_exon.neg.bed.one_line

	awk '{a[$7"\t"$5]++;if (a[$7"\t"$5]==1){if (FNR==1){printf $7"\t"$5"\t"$1"\t"$4"\t"$2"\t"$3}else{print "";printf $7"\t"$5"\t"$1"\t"$4"\t"$2"\t"$3}}else{printf "\t"$2"\t"$3}}END{print ""}' ${sample}.circ2.info.intersect_exon.pos.bed > ${sample}.circ2.info.intersect_exon.pos.bed.one_line



