### Preparatory work ###
#1. Exons in the GTF file
#2. The autosomal genome is converted into one line
#3. N-masked autosomal genome
#4. SNP information of each autosome
		#1	169	170	rs1112347109	G	T
		#1	369	370	rs1108162481	G	A
		#1	481	482	rs1107705765	T	C
		#1	528	529	.	G	A
		#1	554	555	rs1113680398	C	G
		#1	1081	1082	rs1111334233	G	T
#5. Number of autosome

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






	### Sequence alignment
	${bwa_path}/bwa mem -T 19 -a -t 5 \
		${sample}.circRNA_linear_RNA.fa \
		${data_path}/${sample}_1.fq.gz \
		| ${samtools_path}/samtools view -bS - > ${sample}.1.all.bam

	${samtools_path}/samtools view ${sample}.1.all.bam > ${sample}.1.all.sorted.sam.temp_noheader

	${bwa_path}/bwa mem -T 19 -a -t 5 \
		${sample}.circRNA_linear_RNA.fa \
		${data_path}/${sample}_2.fq.gz \
		| ${samtools_path}/samtools view -bS - > ${sample}.2.all.bam

	${samtools_path}/samtools view ${sample}.2.all.bam > ${sample}.2.all.sorted.sam.temp_noheader






	### Perfect alignment reads
	awk '$6!~/[*IDSHP=X]/{print}' ${sample}.1.all.sorted.sam.temp_noheader > ${sample}.1.all.sorted.sam.temp_noheader.complete_map
	awk '$6!~/[*IDSHP=X]/{print}' ${sample}.2.all.sorted.sam.temp_noheader > ${sample}.2.all.sorted.sam.temp_noheader.complete_map

	cut -f1,3,4,6 ${sample}.1.all.sorted.sam.temp_noheader.complete_map | sed "s/M$//g" > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need

	cut -f1,3,4,6 ${sample}.2.all.sorted.sam.temp_noheader.complete_map | sed "s/M$//g" > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need






	### Genome location of read1 which were mapped to linear RNA (forword strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.linear_RNA.bed.pos.one_line ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $(i+1)-$i+1;
			if (a>=$3){
				temp = $3-b+$i-1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = $(i+1) - temp + 1;
				if (m >= $4){
					print temp+$4-1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $(j+1) - $j + 1;
						if (m >= $4){
							print $j"\t"$j+$4-n-1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
				}
				break;
			}
		}
	}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed






	### Genome location of read1 which were mapped to linear RNA (reverse strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.linear_RNA.bed.neg.one_line ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $i-$(i+1)+1;
			if (a>=$3){
				temp = $i-($3-b)+1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = temp - $(i+1) + 1;
				if (m >= $4){
					print temp-$4+1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $j - $(j+1) + 1;
						if (m >= $4){
							print $j"\t"$j-($4-n)+1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
				}
				break;
			}
		}
	}' | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t";for (i=NF;i>5;i--){printf $i"\t"};print $5}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed





	### Genome location of read1 which were mapped to circRNA (forword strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.info.intersect_exon.pos.bed.one_line ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		p=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $(i+1)-$i+1;
			if (a>=$3){
				temp = $3-b+$i-1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = $(i+1) - temp + 1;
				if (m >= $4){
					print temp+$4-1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $(j+1) - $j + 1;
						if (m >= $4){
							print $j"\t"$j+$4-n-1;
							p=1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
					if (p==0){
						for (j=9;j<=NF;j+=2){
							n = m;
							m += $(j+1) - $j + 1;
							if (m >= $4){
								print $j"\t"$j+$4-n-1"\tjunction_reads";
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
				}
				break;
			}
		}

		if (a<$3){
			raw_3 = $3;
			$3 = $3 - a;
			a=0;
			for (i=9;i<=NF;i+=2){
				b = a;
				a += $(i+1)-$i+1;
				if (a>=$3){
					temp = $3-b+$i-1;
					printf $1"\t"$2"\t"raw_3"\t"$4"\t"temp"\t";
					n = m;
					m = $(i+1) - temp + 1;
					if (m >= $4){
						print temp+$4-1;
					}else{
						printf $(i+1)"\t";
						for (j=i+2;j<=NF;j+=2){
							n = m;
							m += $(j+1) - $j + 1;
							if (m >= $4){
								print $j"\t"$j+$4-n-1;
								p=1;
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
					break;
				}
			}
		}
		
	}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed





	### Genome location of read1 which were mapped to circRNA (reverse strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.info.intersect_exon.neg.bed.one_line ${sample}.1.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		p=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $i-$(i+1)+1;
			if (a>=$3){
				temp = $i-($3-b)+1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = temp - $(i+1) + 1;
				if (m >= $4){
					print temp-$4+1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $j - $(j+1) + 1;
						if (m >= $4){
							print $j"\t"$j-($4-n)+1;
							p=1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
					
					if (p==0){
						for (j=9;j<=NF;j+=2){
							n = m;
							m += $j - $(j+1) + 1;
							if (m >= $4){
								print $j"\t"$j-($4-n)+1"\tjunction_reads";
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}				
					
				}
				break;
			}
		}
			
		if (a<$3){
			raw_3 = $3;
			$3 = $3 - a;
			a=0;
			for (i=9;i<=NF;i+=2){
				b = a;
				a += $i-$(i+1)+1;
				if (a>=$3){
					temp = $i-($3-b)+1;
					printf $1"\t"$2"\t"raw_3"\t"$4"\t"temp"\t";
					n = m;
					m = temp - $(i+1) + 1;
					if (m >= $4){
						print temp-$4+1;
					}else{
						printf $(i+1)"\t";
						for (j=i+2;j<=NF;j+=2){
							n = m;
							m += $j - $(j+1) + 1;
							if (m >= $4){
								print $j"\t"$j-($4-n)+1;
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
					break;
				}
			}
		}
	}' | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t";if (NF%2==0){for (i=NF;i>5;i--){printf $i"\t"};print $5}else{for (i=NF-1;i>=5;i--){printf $i"\t"};print $NF}}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed





	### Genome location of read2 which were mapped to linear RNA (forword strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.linear_RNA.bed.pos.one_line ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $(i+1)-$i+1;
			if (a>=$3){
				temp = $3-b+$i-1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = $(i+1) - temp + 1;
				if (m >= $4){
					print temp+$4-1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $(j+1) - $j + 1;
						if (m >= $4){
							print $j"\t"$j+$4-n-1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
				}
				break;
			}
		}
	}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed





	### Genome location of read2 which were mapped to linear RNA (reverse strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.linear_RNA.bed.neg.one_line ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $i-$(i+1)+1;
			if (a>=$3){
				temp = $i-($3-b)+1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = temp - $(i+1) + 1;
				if (m >= $4){
					print temp-$4+1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $j - $(j+1) + 1;
						if (m >= $4){
							print $j"\t"$j-($4-n)+1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
				}
				break;
			}
		}
	}' | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t";for (i=NF;i>5;i--){printf $i"\t"};print $5}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed





	### Genome location of read2 which were mapped to circRNA (forword strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.info.intersect_exon.pos.bed.one_line ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		p=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $(i+1)-$i+1;
			if (a>=$3){
				temp = $3-b+$i-1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = $(i+1) - temp + 1;
				if (m >= $4){
					print temp+$4-1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $(j+1) - $j + 1;
						if (m >= $4){
							print $j"\t"$j+$4-n-1;
							p=1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
					if (p==0){
						for (j=9;j<=NF;j+=2){
							n = m;
							m += $(j+1) - $j + 1;
							if (m >= $4){
								print $j"\t"$j+$4-n-1"\tjunction_reads";
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
				}
				break;
			}
		}

		if (a<$3){
			raw_3 = $3;
			$3 = $3 - a;
			a=0;
			for (i=9;i<=NF;i+=2){
				b = a;
				a += $(i+1)-$i+1;
				if (a>=$3){
					temp = $3-b+$i-1;
					printf $1"\t"$2"\t"raw_3"\t"$4"\t"temp"\t";
					n = m;
					m = $(i+1) - temp + 1;
					if (m >= $4){
						print temp+$4-1;
					}else{
						printf $(i+1)"\t";
						for (j=i+2;j<=NF;j+=2){
							n = m;
							m += $(j+1) - $j + 1;
							if (m >= $4){
								print $j"\t"$j+$4-n-1;
								p=1;
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
					break;
				}
			}
		}
		
	}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed





	### Genome location of read2 which were mapped to circRNA (reverse strand)
	awk 'NR==FNR{a[$2]=$0}NR>FNR&&a[$2]{print $0"\t"a[$2]}' ${sample}.circ2.info.intersect_exon.neg.bed.one_line ${sample}.2.all.sorted.sam.temp_noheader.complete_map.need | awk '{
		a=0;
		m=0;
		p=0;
		for (i=9;i<=NF;i+=2){
			b = a;
			a += $i-$(i+1)+1;
			if (a>=$3){
				temp = $i-($3-b)+1;
				printf $1"\t"$2"\t"$3"\t"$4"\t"temp"\t";
				n = m;
				m = temp - $(i+1) + 1;
				if (m >= $4){
					print temp-$4+1;
				}else{
					printf $(i+1)"\t";
					for (j=i+2;j<=NF;j+=2){
						n = m;
						m += $j - $(j+1) + 1;
						if (m >= $4){
							print $j"\t"$j-($4-n)+1;
							p=1;
							break;
						}else{
							printf $j"\t"$(j+1)"\t";
						}
					}
					
					if (p==0){
						for (j=9;j<=NF;j+=2){
							n = m;
							m += $j - $(j+1) + 1;
							if (m >= $4){
								print $j"\t"$j-($4-n)+1"\tjunction_reads";
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}				
					
				}
				break;
			}
		}
			
		if (a<$3){
			raw_3 = $3;
			$3 = $3 - a;
			a=0;
			for (i=9;i<=NF;i+=2){
				b = a;
				a += $i-$(i+1)+1;
				if (a>=$3){
					temp = $i-($3-b)+1;
					printf $1"\t"$2"\t"raw_3"\t"$4"\t"temp"\t";
					n = m;
					m = temp - $(i+1) + 1;
					if (m >= $4){
						print temp-$4+1;
					}else{
						printf $(i+1)"\t";
						for (j=i+2;j<=NF;j+=2){
							n = m;
							m += $j - $(j+1) + 1;
							if (m >= $4){
								print $j"\t"$j-($4-n)+1;
								break;
							}else{
								printf $j"\t"$(j+1)"\t";
							}
						}
					}
					break;
				}
			}
		}
	}' | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t";if (NF%2==0){for (i=NF;i>5;i--){printf $i"\t"};print $5}else{for (i=NF-1;i>=5;i--){printf $i"\t"};print $NF}}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed





	### Combine genome location of forword and reverse strand
	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed

	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed

	cat ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed

	cat ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed

	rm ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.neg.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.pos.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.neg.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.pos.genome.bed





	### read1 and read2 align to the same transcript simultaneously
	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$1"\t"$2]{print}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair

	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$1"\t"$2]{print}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair

	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$1"\t"$2]{print}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair

	awk 'NR==FNR{a[$1"\t"$2]++}NR>FNR&&a[$1"\t"$2]{print}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair

	rm ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed






	### The location of read1 and read2 alignment to the transcriptome maps to a unique location on the genome
	cut -f1,5- ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair | awk '{a[$0]++}END{for (i in a)print i}' | cut -f1 | awk '{a[$0]++}END{for (i in a){if (a[i]>1){print i}}}' > ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq

	cut -f1,5- ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair | awk '{a[$0]++}END{for (i in a)print i}' | cut -f1 | awk '{a[$0]++}END{for (i in a){if (a[i]>1){print i}}}' >> ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq


	cut -f1,5- ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair | awk '{a[$0]++}END{for (i in a)print i}' | cut -f1 | awk '{a[$0]++}END{for (i in a){if (a[i]>1){print i}}}' > ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq

	cut -f1,5- ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair | awk '{a[$0]++}END{for (i in a)print i}' | cut -f1 | awk '{a[$0]++}END{for (i in a){if (a[i]>1){print i}}}' >> ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq

	if [[ -s ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq ]]
	then

		awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci

		awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci

	else

		cp ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci

		cp ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci

	fi


	if [[ -s ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq ]]
	then

		awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci

		awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1]{print $0}'  ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci

	else

		cp ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci

		cp ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci

	fi

	rm ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.not_uniq ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.not_uniq ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair

	sed "s/_/\t/g" ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci | cut -f1,2,6- | awk '{a[$0]++}END{for (i in a)print i}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci

	sed "s/_/\t/g" ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci | cut -f1,2,6- | awk '{a[$0]++}END{for (i in a)print i}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci

	sed "s/_/\t/g" ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci | cut -f1,2,6- | awk '{a[$0]++}END{for (i in a)print i}' > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci

	sed "s/_/\t/g" ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci | cut -f1,2,6- | awk '{a[$0]++}END{for (i in a)print i}' > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci


	rm ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci





	### The reads are not matched to circular RNA and linear RNA simultaneously
	awk 'NR==FNR{a[$1]++}NR>FNR&&a[$1]{print $1}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci > ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA

	awk 'NR==FNR{a[$1]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear

	awk 'NR==FNR{a[$1]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear

	awk 'NR==FNR{a[$1]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear

	awk 'NR==FNR{a[$1]++}NR>FNR&&!a[$1]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear

	rm ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.in_linear_RNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci





	### reads aligned to linear RNA 
	cut -f5,8 ${sample}.circ2.info.intersect_exon.bed | sed "s/_[0-9]*//g" | sort | uniq | sed "s/:/@/g;s/|/@/g" > circRNA_gene_pair

	awk '{a[$2]++;if (a[$2]==1){b[$2]=$1}else{b[$2]=b[$2]"_______"$1}}END{for (i in b)print i"\t"b[i]}' circRNA_gene_pair | awk 'NR==FNR{a[$1]=$2}NR>FNR&&a[$2]{print $0"\t"a[$2]}' - ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear | awk -F "\t|_______" '{for (i=5;i<=NF;i++){if ($i~/@/){p=i;break}};for (i=p;i<=NF;i++){for (j=1;j<p;j++){printf $j"\t"};print $i}}' | sed "s/@/\t/g" > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA

	awk '{a[$2]++;if (a[$2]==1){b[$2]=$1}else{b[$2]=b[$2]"_______"$1}}END{for (i in b)print i"\t"b[i]}' circRNA_gene_pair | awk 'NR==FNR{a[$1]=$2}NR>FNR&&a[$2]{print $0"\t"a[$2]}' - ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear | awk -F "\t|_______" '{for (i=5;i<=NF;i++){if ($i~/@/){p=i;break}};for (i=p;i<=NF;i++){for (j=1;j<p;j++){printf $j"\t"};print $i}}' | sed "s/@/\t/g" > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA





	### Position relationship between circular RNA and reads
	awk '{if ($(NF-3)<$(NF-1)){print $1"\t"$2"\tleft\t"$(NF-2)":"$(NF-1)"|"$NF}else if ($3>$NF){print $1"\t"$2"\tright\t"$(NF-2)":"$(NF-1)"|"$NF}else{print $0}}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right

	awk '{if ($(NF-3)<$(NF-1)){print $1"\t"$2"\tleft\t"$(NF-2)":"$(NF-1)"|"$NF}else if ($3>$NF){print $1"\t"$2"\tright\t"$(NF-2)":"$(NF-1)"|"$NF}else{print $0}}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right

	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right | grep -E "left|right" | awk '{a[$0]++}END{for (i in a){if (a[i]>1){print i}}}' | awk '{print $1"\t"$2"\t"$NF}' > ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.both_out_circRNA





	### reads that cross or cover circular RNA
	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.both_out_circRNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross

	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.both_out_circRNA ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross

	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross | grep -E "left|right" | awk '{print $1"\t"$2"\t"$NF}' | sort | uniq -d > ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cross_circRNA





	### reads that cover circular RNA
	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cross_circRNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover

	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cross_circRNA ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover





	### read1 covers the region of circular RNA
	awk '{
		temp=0;
		read_length=0;
		if ($3=="left"||$3=="right"){
			print $0;
		}else{
			m=0;
			n=0;
			x=0;
			y=0;
			for (i=3;i<=NF-3;i+=2){
				read_length+=$(i+1)-$i+1;
				if ($(NF-1)>=$i && $(i+1)>=$(NF-1)){
					m=i;
				}
				if ($NF>=$i && $(i+1)>=$NF){
					n=i;
				}
			}
			for (i=4;i<=NF-4;i+=2){
				if ($(NF-1)>$i && $(i+1)>$(NF-1)){
					x=i+1;
				}
				if ($NF>$i && $(i+1)>$NF){
					y=i+1;
				}
			}
			if (x==0 && y==0){
				if (m==0){
					if (n==0){
						printf $1"\t"$2"\t"read_length"\t1\t";
						for (i=3;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;	
					}else{
						printf $1"\t"$2"\t"read_length"\t1\t";
						for (i=3;i<n;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;			
					}		
				}else{
					if (n==0){
						for (i=3;i<m;i+=2){
							temp+=$(i+1)-$i+1;
						}
						temp+=$(NF-1)-$m+1;
						printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
						for (i=m+2;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;
					}else{
						for (i=3;i<m;i+=2){
							temp+=$(i+1)-$i+1;
						}
						temp+=$(NF-1)-$m+1;
						if (m==n){
							print $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
						}else{
							printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
							for (i=m+2;i<n;i+=2){
								printf $i"\t"$(i+1)"\t";
							}
							print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
						}
					}
				}
			}else if (x==0 && y!=0){
				if (m==0){
					printf $1"\t"$2"\t"read_length"\t1\t";
					for (i=3;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;				
				}else{
					for (i=3;i<m;i+=2){
						temp+=$(i+1)-$i+1;
					}
					temp+=$(NF-1)-$m+1;
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
					for (i=m+2;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;			
				}
			}else if (x!=0 && y==0){
				if (n==0){
						for (i=3;i<x;i+=2){
							temp+=$(i+1)-$i+1;
						}
						printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
						for (i=x+2;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;			
				}else{
					for (i=3;i<x;i+=2){
						temp+=$(i+1)-$i+1;
					}
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
					for (i=x+2;i<n;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					if (x==n){
						print $(NF-2)":"$(NF-1)"|"$NF;
					}else{
						print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
					}
				}
			}else{
				if (x==y){
					print $0"\tcross";
				}else{	
					for (i=3;i<x;i+=2){
						temp+=$(i+1)-$i+1;
					}
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
					for (i=x+2;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;
				}
			}
		}
	}' ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region






	### read2 covers the region of circular RNA
	awk '{
		temp=0;
		read_length=0;
		if ($3=="left"||$3=="right"){
			print $0;
		}else{
			m=0;
			n=0;
			x=0;
			y=0;
			for (i=3;i<=NF-3;i+=2){
				read_length+=$(i+1)-$i+1;
				if ($(NF-1)>=$i && $(i+1)>=$(NF-1)){
					m=i;
				}
				if ($NF>=$i && $(i+1)>=$NF){
					n=i;
				}
			}
			for (i=4;i<=NF-4;i+=2){
				if ($(NF-1)>$i && $(i+1)>$(NF-1)){
					x=i+1;
				}
				if ($NF>$i && $(i+1)>$NF){
					y=i+1;
				}
			}
			if (x==0 && y==0){
				if (m==0){
					if (n==0){
						printf $1"\t"$2"\t"read_length"\t1\t";
						for (i=3;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;	
					}else{
						printf $1"\t"$2"\t"read_length"\t1\t";
						for (i=3;i<n;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;			
					}		
				}else{
					if (n==0){
						for (i=3;i<m;i+=2){
							temp+=$(i+1)-$i+1;
						}
						temp+=$(NF-1)-$m+1;
						printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
						for (i=m+2;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;
					}else{
						for (i=3;i<m;i+=2){
							temp+=$(i+1)-$i+1;
						}
						temp+=$(NF-1)-$m+1;
						if (m==n){
							print $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
						}else{
							printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
							for (i=m+2;i<n;i+=2){
								printf $i"\t"$(i+1)"\t";
							}
							print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
						}
					}
				}
			}else if (x==0 && y!=0){
				if (m==0){
					printf $1"\t"$2"\t"read_length"\t1\t";
					for (i=3;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;				
				}else{
					for (i=3;i<m;i+=2){
						temp+=$(i+1)-$i+1;
					}
					temp+=$(NF-1)-$m+1;
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$(NF-1)"\t"$(m+1)"\t";
					for (i=m+2;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;			
				}
			}else if (x!=0 && y==0){
				if (n==0){
						for (i=3;i<x;i+=2){
							temp+=$(i+1)-$i+1;
						}
						printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
						for (i=x+2;i<=NF-3;i+=2){
							printf $i"\t"$(i+1)"\t";
						}
						print $(NF-2)":"$(NF-1)"|"$NF;			
				}else{
					for (i=3;i<x;i+=2){
						temp+=$(i+1)-$i+1;
					}
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
					for (i=x+2;i<n;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					if (x==n){
						print $(NF-2)":"$(NF-1)"|"$NF;
					}else{
						print $n"\t"$NF"\t"$(NF-2)":"$(NF-1)"|"$NF;
					}
				}
			}else{
				if (x==y){
					print $0"\tcross";
				}else{	
					for (i=3;i<x;i+=2){
						temp+=$(i+1)-$i+1;
					}
					printf $1"\t"$2"\t"read_length"\t"temp"\t"$x"\t"$(x+1)"\t";
					for (i=x+2;i<y;i+=2){
						printf $i"\t"$(i+1)"\t";
					}
					print $(NF-2)":"$(NF-1)"|"$NF;
				}
			}
		}
	}' ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region

	cat ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region | awk '$NF=="cross"{print $1"\t"$2}' | sort |uniq > ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.cross





	### At least 10bp across the back splicing site
	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2]' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.cross ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region | grep -v -E "left|right" | awk '{a=0;for (i=5;i<=NF-1;i+=2){a+=$(i+1)-$i+1};print $1"\t"$2"\t"$4-1"\t"a"\t"$3-$4+1-a"\t"$NF}' | awk '($3>=10||$5>=10)&&$4>=10{print $0}' > ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.valid

	awk 'NR==FNR{a[$0]++}NR>FNR&&!a[$1"\t"$2]' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.cross ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region | grep -v -E "left|right" | awk '{a=0;for (i=5;i<=NF-1;i+=2){a+=$(i+1)-$i+1};print $1"\t"$2"\t"$4-1"\t"a"\t"$3-$4+1-a"\t"$NF}' | awk '($3>=10||$5>=10)&&$4>=10{print $0}' >> ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.valid

	awk 'NR==FNR{a[$1"\t"$2"\t"$NF]++}NR>FNR&&a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.valid ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region

	awk 'NR==FNR{a[$1"\t"$2"\t"$NF]++}NR>FNR&&a[$1"\t"$2"\t"$NF]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.valid ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.valid.region

	rm ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.left_midlle_right ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.both_out_circRNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover_or_cross ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cross_circRNA ${sample}.1.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region ${sample}.2.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.cross ${sample}.all.sorted.sam.temp_noheader.complete_map.linear_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.circRNA.cover.region.valid





	### Reads on circular RNA
	cat ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear | awk '{a=0;b=0;p=0;for (i=3;i<=NF-4;i+=2){a+=$(i+1)-$i+1;if ($i>$(i+2)){p=i+2;break}};for (i=p;i<=NF-2;i+=2){b+=$(i+1)-$i+1};print $1"\t"$2"\t"a"\t"b}' | awk '$3>=10&&$4>=10{print}' | cut -f1,2 | sort |uniq > ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid

	awk 'NR==FNR{a[$0]++}NR>FNR&&a[$1"\t"$2]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear > ${sample}.1.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid

	awk 'NR==FNR{a[$0]++}NR>FNR&&a[$1"\t"$2]{print $0}' ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear > ${sample}.2.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid

	rm ${sample}.all.sorted.sam.temp_noheader.complete_map.circular_RNA.genome.bed.pair.only_loci.loci.uniq_in_circ_linear.valid






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

