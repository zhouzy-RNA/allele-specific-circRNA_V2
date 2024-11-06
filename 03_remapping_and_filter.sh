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
#1. Read1 and read2 were mapped to the reference transcriptome, respectively 
#2. Perfect alignment reads
#3. The transcriptome location of reads were mapped to the genome location.
#4. filter
#	a. read1 and read2 align to the same transcript simultaneously
#	b. The genome location of read1 and read2 were uniq
#	c. The reads are not matched to circular RNA and linear RNA simultaneously
#	d. At least 10bp across the back splicing site





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


