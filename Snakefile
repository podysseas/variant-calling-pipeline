# workflow implementation with Snakemake 

# pwd : ~/dc_workshop/my_workshop/
# mkdir <files> : data , results/bam , resuls/bcf , results/sam , results,bam
# the above files should be empty
# ecoli should be in the <pwd> file 

rule all:
 	input :
 		"data/ecoli_rel606.fasta",
 		"trimmed_fastq/SRR2584866_1.trim.fastq.gz",
		"trimmed_fastq/SRR2584866_1un.trim.fastq.gz",
		"trimmed_fastq/SRR2584866_2.trim.fastq.gz",
		"trimmed_fastq/SRR2584866_2un.trim.fastq.gz",
		"trimmed_fastq/SRR2584863_1.trim.fastq.gz",
		"trimmed_fastq/SRR2584863_1un.trim.fastq.gz",
		"trimmed_fastq/SRR2584863_2.trim.fastq.gz",
		"trimmed_fastq/SRR2584863_2un.trim.fastq.gz",
		"trimmed_fastq/SRR2589044_1.trim.fastq.gz",
		"trimmed_fastq/SRR2589044_1un.trim.fastq.gz",
		"trimmed_fastq/SRR2589044_2.trim.fastq.gz",
		"trimmed_fastq/SRR2589044_2un.trim.fastq.gz",
		"trimmed_fastq/SRR2584866_1.trim.fastq",
		"trimmed_fastq/SRR2584866_1un.trim.fastq",
		"trimmed_fastq/SRR2584866_2.trim.fastq",
		"trimmed_fastq/SRR2584866_2un.trim.fastq",
		"trimmed_fastq/SRR2584863_1.trim.fastq",
		"trimmed_fastq/SRR2584863_1un.trim.fastq",
		"trimmed_fastq/SRR2584863_2.trim.fastq",
		"trimmed_fastq/SRR2584863_2un.trim.fastq",
		"trimmed_fastq/SRR2589044_1.trim.fastq",
		"trimmed_fastq/SRR2589044_1un.trim.fastq",
		"trimmed_fastq/SRR2589044_2.trim.fastq",
		"trimmed_fastq/SRR2589044_2un.trim.fastq",
		# "results/bcf/SRR2589044_raw.bcf",
		# "results/bcf/SRR2584866_raw.bcf",
		# "results/bcf/SRR2584863_raw.bcf",
		# "results/sam/SRR2584863.aligned.sam",
		# "results/sam/SRR2584866.aligned.sam",
		# "results/sam/SRR2589044.aligned.sam",
		# "results/bam/SRR2584866.aligned.bam",
		# "results/bam/SRR2584863.aligned.bam",
		# "results/bam/SRR2589044.aligned.bam",     
		# "results/bam/SRR2584866.aligned.sorted.bam",
		# "results/bam/SRR2584863.aligned.sorted.bam",
		# "results/bam/SRR2589044.aligned.sorted.bam",
		"results/vcf/SRR2589044_final_variants.vcf",
		"results/vcf/SRR2584863_final_variants.vcf",
		"results/vcf/SRR2584866_final_variants.vcf"
rule trimming:
	input:
		first="untrimmed_fastq/{sample}_1.fastq.gz"
		#second="untrimmed_fastq/{sample}_2.fastq.gz"
	output:
		first_trimmed="trimmed_fastq/{sample}_1.trim.fastq.gz",
		first_untrimmed="trimmed_fastq/{sample}_1un.trim.fastq.gz",
		second_trimmed="trimmed_fastq/{sample}_2.trim.fastq.gz",
		second_untrimmed="trimmed_fastq/{sample}_2un.trim.fastq.gz"
	priority:40
	shell: 
		"trimmomatic PE {input.first} untrimmed_fastq/{wildcards.sample}_2.fastq.gz \
		{output.first_trimmed} {output.first_untrimmed} {output.second_trimmed} {output.second_untrimmed}\
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15"
    #  use a sliding window of size 4 that will remove bases if their phred 
    #  score is below 20 and discard any reads that do not have at least 25 bases 
    #  remaining after this trimming step.
rule indexing:                                     
	input: 
		gen="ecoli_rel606.fasta"
	output:
		indexed="data/ecoli_rel606.fasta" # problem with cyclic dependency,that's why should be in another folder 
	priority: 50
	shell:
		"bwa index {input.gen} > {output.indexed}"
rule gun_zipping:								   
	input: 
		file="trimmed_fastq/{G_SAMPLES}.trim.fastq.gz" # here there is an error not finding one file 
		# not problem with implementation ,maybe lacks in speed in case of using multiple threads
	output:
		"trimmed_fastq/{G_SAMPLES}.trim.fastq"
	priority:30
	shell:
		"gunzip {input.file}"

rule align_reads:
	input: 
		first="trimmed_fastq/{SAMPLES}_1.trim.fastq",
		gen="ecoli_rel606.fasta"
	output:
		"results/sam/{SAMPLES}.aligned.sam"
	priority:20
	shell:
		 "bwa mem {input.gen} {input.first} trimmed_fastq/{wildcards.SAMPLES}_2.trim.fastq > {output}"
rule convert_to_sam:
	input:
		"results/sam/{C_SAMPLES}.aligned.sam"
	output:
		"results/bam/{C_SAMPLES}.aligned.bam"
	priority:10
	shell:
		"samtools view -S -b {input} > {output}"
rule sorting:
	input:
		"results/bam/{S_SAMPLES}.aligned.bam"
	output:
		"results/bam/{S_SAMPLES}.aligned.sorted.bam"
	priority:5
	shell:
		"samtools sort -o {output} {input}"
rule variant_calling_first_step:
	input:
		gen="ecoli_rel606.fasta",
		sort="results/bam/{V_SAMPLES}.aligned.sorted.bam"
	output:
		bcf="results/bcf/{V_SAMPLES}_raw.bcf"
	priority:4
	shell:
		"bcftools mpileup -O b -o {output.bcf} -f {input.gen} {input.sort}"
rule detect_polymorphisms:
	input:
		bcf="results/bcf/{V_SAMPLES}_raw.bcf"
	output:
		variants="results/bcf/{V_SAMPLES}_variants.vcf"
	priority:3
	shell:
		"bcftools call --ploidy 1 -m -v -o {output.variants} {input.bcf}"
rule Filter_and_report:
	input:
		variants="results/bcf/{V_SAMPLES}_variants.vcf"
	output:
		f_variants="results/vcf/{V_SAMPLES}_final_variants.vcf"
	priority:2
	shell:										
		"vcfutils.pl varFilter {input.variants}  > {output.f_variants}"
