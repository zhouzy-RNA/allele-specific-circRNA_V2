chr_number=18
chr_seq="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sscrofa11.1.dna.chr.fa"
gtf="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sus_scrofa.Sscrofa11.1.90.chromsome.gtf"
vcf_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/vcf"
script_path="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/script"
N_masked_ref_seq="/data_group/zhouzhongyin/zhouzhongyin1/circRNA/allele_specific_circRNA/pig/ref/Sscrofa11.1.dna.chr.masked.fa"
bwa_path="/opt/service/fermi.kit"

### Build index
${bwa_path}/bwa index -a bwtsw ${chr_seq}



### The chromosome sequence is arranged in a row
perl ${script_path}/FastaGenomeToNewFormat.pl \
	-origin_fasta ${chr_seq} \
	-seq1row_fasta ${chr_seq}.1row



### Split chromosomes and generate N-masked chromosome sequences
for((i=1;i<=${chr_num};i++));
do

	sed -n "$[i*2-1],$[i*2]p" ${chr_seq}.fa.1row > ${chr_seq}.chr${line}.fa.1row;

	perl ${script_path}/GenomeSequenceMaskN.pl \
		-input_genome ${chr_seq}.chr${line}.fa.1row \
		-input_snp ${vcf_path}/chr${line}_overlapping.chr_posit.vcf \
		-output_maskedgenome ${chr_seq}.chr${line}.masked.fa.1row

	cat ${chr_seq}.chr${line}.masked.fa.1row >> ${N_masked_ref_seq}

done


### Exons in the GTF file
grep $'\t'exon$'\t' $gtf | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$18}' | sed "s/\"//g;s/;//g" > ${gtf}.exon.bed





