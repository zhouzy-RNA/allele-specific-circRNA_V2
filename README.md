# allele-specific-circRNA_V2
This pipline is used to identify allele-specific circRNAs within individual.

1. The preparation required to identify allele-specific circRNAs is in 01_Preparatory_work.sh.
2. Allele_specific_circRNA.all_process.sh contains 02_circRNA_annotation.sh, 03_remapping_and_filter.sh and 04_allele_specific_circRNA.sh.
3. Using this pipeline only requires modifying the the path of data, software, and the results.
4. Result description: 
	1 column: Sample Name
	2 column: circRNA
	3 column: SNP Position
	4 column: SNP ID
	5 column: Ref base
	6 column: Alt base
	7 column: n1
	8 column: n2
	9 column: n3
	10 column: n4
	11 column: Odd Ratio
	12 column: lower of 95% confidence interval
	13 column: upper of 95% confidence
