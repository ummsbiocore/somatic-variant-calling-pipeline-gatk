$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			def recursive = (type == "folder") ? "--recursive" : ""
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.genome){params.genome = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.known_snps){params.known_snps = ""} 
if (!params.known_indels){params.known_indels = ""} 
if (!params.known_snps_index){params.known_snps_index = ""} 
if (!params.known_indels_index){params.known_indels_index = ""} 
if (!params.groups_file){params.groups_file = ""} 
if (!params.compare_file){params.compare_file = ""} 
if (!params.annotation){params.annotation = ""} 
if (!params.annotation_index){params.annotation_index = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)

g_2_0_g_0 = file(params.genome, type: 'any')
g_2_0_g_10 = file(params.genome, type: 'any')
g_2_1_g_23 = file(params.genome, type: 'any')
if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_4_0_g_3}
  } else {  
	g_4_0_g_3 = Channel.empty()
 }

Channel.value(params.mate).set{g_5_1_g_3}
g_11_2_g_9 = file(params.known_snps, type: 'any')
g_12_3_g_9 = file(params.known_indels, type: 'any')
g_15_4_g_9 = file(params.known_snps_index, type: 'any')
g_16_5_g_9 = file(params.known_indels_index, type: 'any')
g_18_1_g_17 = file(params.groups_file, type: 'any')
g_19_2_g_17 = file(params.compare_file, type: 'any')
g_24_2_g_23 = file(params.annotation, type: 'any')
g_27_3_g_23 = file(params.annotation_index, type: 'any')

build_BWA_index = params.Check_Build_BWA.build_BWA_index
//* params.bwa_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process Check_Build_BWA {

input:
 path genome

output:
 path "$index"  ,emit:g_0_bwaindex00_g_1 

stageInMode 'copy'

when:
build_BWA_index == true && ((params.run_BWA && (params.run_BWA == "yes")) || !params.run_BWA)

script:

bwa_build_parameters = params.Check_Build_BWA.bwa_build_parameters
basename = genome.baseName
index = "BWAIndex" 


"""
mkdir -p $index && mv $genome $index/. && cd $index
bwa index ${bwa_build_parameters} ${genome} 
samtools faidx ${genome} 
gatk CreateSequenceDictionary -R ${genome} 
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process check_BWA_files {

input:
 path bwa

output:
 path "*/${bwa2}" ,optional:true  ,emit:g_1_bwaindex02_g_3 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
stageInMode 'copy'

when:
(params.run_BWA && (params.run_BWA == "yes")) || !params.run_BWA

script:
(cmd, bwa2) = pathChecker(bwa, params.bwa_index, "folder")
"""
$cmd
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 8
    $MEMORY = 32
}
//* platform
//* platform
//* autofill

process bwa_align {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_aligned_reads.bam$/) "bwa/$filename"}
input:
 tuple val(name),file(reads)
 val mate
 path ref

output:
 tuple val(name), file("${name}_aligned_reads.bam")  ,emit:g_3_bam_file00_g_6 

script:

readGroup = "@RG\\tID:${name}\\tLB:${name}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${name}"
bwa_mem_align_parameters = params.bwa_align.bwa_mem_align_parameters

"""
indexname=\$(ls ${ref}/*.fa )
ls $ref
bwa mem -t ${task.cpus} ${bwa_mem_align_parameters} \
	-R '${readGroup}' \
	-P \${indexname}  \
	${reads} \
	| samtools view -bh - | samtools sort - -o ${name}_aligned_reads.bam
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process markDuplicates {

input:
 tuple val(name), file(reads)

output:
 tuple val(name), file("${name}_sorted_dedup.bam")  ,emit:g_6_mapped_reads00_g_9 
 path "${name}_dedup_metrics.txt"  ,emit:g_6_txtFile11 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
mkdir -p tmp/${name}
gatk --java-options '-Djava.io.tmpdir=tmp/${name}' \
 MarkDuplicates \
-I ${reads} \
-M ${name}_dedup_metrics.txt \
-O ${name}_sorted_dedup.bam 
rm -r tmp/${name}
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process build_gatk4_genome_dictionary {

input:
 path genome

output:
 path "genomeDict"  ,emit:g_10_genomeDict01_g_9 

container 'quay.io/ummsbiocore/gatk:1.0.0'
stageInMode 'copy'

script:

basename = genome.baseName

"""
mkdir -p genomeDict && mv ${genome} genomeDict/ && cd genomeDict
samtools faidx ${genome}
gatk CreateSequenceDictionary -R ${genome} 
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 2
}
//* platform
//* platform
//* autofill

process Base_Recalibration {

input:
 tuple val(name), file(input_bam)
 path reference
 path known_snps
 path known_indels
 path known_snps_index
 path known_indels_index

output:
 tuple val("${name}"), file("${name}_recal_data.txt")  ,emit:g_9_outputFileTxt02_g_13 
 tuple val("${name}"), file("${input_bam}")  ,emit:g_9_bamFile10_g_13 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
basecal_parameters = params.Base_Recalibration.basecal_parameters

"""
indexname=\$(ls ${reference}/*.fa )

gatk BaseRecalibrator \
	-R \${indexname}  ${basecal_parameters} \
	-I ${input_bam} \
	--known-sites ${known_snps} \
    --known-sites ${known_indels} \
	-O ${name}_recal_data.txt
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 5
}
//* platform
//* platform
//* autofill

process BQSR {

input:
 tuple val(name), file(bam)
 path reference
 tuple val(name), file(recalibration_table)

output:
 path "${name}_recal.bam*"  ,emit:g_13_bamFile00_g_17 

container 'quay.io/ummsbiocore/gatk:1.0.0'

script:

"""
indexname=\$(ls ${reference}/*.fa )
gatk ApplyBQSR \
    -R \$indexname \
    -I ${bam} \
    -bqsr ${recalibration_table} \
    -O ${name}_recal.bam
    
samtools index ${name}_recal.bam
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process prepare_comparisons {

input:
 path bam
 path groups
 path comparisons

output:
 path "*"  ,emit:g_17_inputDir00_g_21 

container 'python'

script:
"""
prepare_comparisons.py --groups ${groups} --comparisons ${comparisons}
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
//* platform
//* autofill

process mutect2 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${comparison}.vcf.gz$/) "mutect2/$filename"}
input:
 path comparison
 path reference

output:
 path "${comparison}.vcf.gz"  ,emit:g_21_vcfFile00_g_23 

container 'quay.io/ummsbiocore/gatk:1.0.0'

script:

"""
indexname=\$(ls ${reference}/*.fa )

controls=\$(ls ${comparison}/controls/*.bam)
control_string=\$""

for i in \$controls
do
    control_string=\$control_string" -I \$i"
done

treatments=\$(ls ${comparison}/treatments/*.bam)
treatment_string=\$""

for i in \$treatments
do
    treatment_string=\$treatment_string" -I \$i"
done

normal_string=\$""

for i in \$controls
do
    normal_string=\$normal_string" --normal-sample \$(basename \$i _recal.bam)"
done

gatk Mutect2 \
	-R \$indexname \
	\$control_string \
	\$treatment_string \
	\$normal_string \
	-O ${comparison}.vcf.gz
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process vep {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${basename}.*$/) "vep/$filename"}
input:
 path input_vcf
 path genome
 path annotation
 path annotation_index

output:
 path "${basename}*"  ,emit:g_23_vcfFile00 

container 'quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1'
disk 1000.GB

script:
basename = input_vcf.baseName.replaceAll(".vcf", "")
annotation_name = annotation.baseName

"""
vep --fasta ${genome} --gtf ${annotation} -i ${input_vcf} -o ${basename} --symbol --biotype
sed -i 's/http://g' ${basename}_summary.html

mv ${basename} ${basename}.vcf
"""
}


workflow {


Check_Build_BWA(g_2_0_g_0)
g_0_bwaindex00_g_1 = Check_Build_BWA.out.g_0_bwaindex00_g_1

g_0_bwaindex00_g_1= g_0_bwaindex00_g_1.ifEmpty(ch_empty_file_1) 


if (!((params.run_BWA && (params.run_BWA == "yes")) || !params.run_BWA)){
g_0_bwaindex00_g_1.set{g_1_bwaindex02_g_3}
} else {

check_BWA_files(g_0_bwaindex00_g_1)
g_1_bwaindex02_g_3 = check_BWA_files.out.g_1_bwaindex02_g_3
}


bwa_align(g_4_0_g_3,g_5_1_g_3,g_1_bwaindex02_g_3)
g_3_bam_file00_g_6 = bwa_align.out.g_3_bam_file00_g_6


markDuplicates(g_3_bam_file00_g_6)
g_6_mapped_reads00_g_9 = markDuplicates.out.g_6_mapped_reads00_g_9
g_6_txtFile11 = markDuplicates.out.g_6_txtFile11


build_gatk4_genome_dictionary(g_2_0_g_10)
g_10_genomeDict01_g_9 = build_gatk4_genome_dictionary.out.g_10_genomeDict01_g_9
(g_10_genomeDict01_g_13,g_10_genomeDict01_g_21) = [g_10_genomeDict01_g_9,g_10_genomeDict01_g_9]


Base_Recalibration(g_6_mapped_reads00_g_9,g_10_genomeDict01_g_9,g_11_2_g_9,g_12_3_g_9,g_15_4_g_9,g_16_5_g_9)
g_9_outputFileTxt02_g_13 = Base_Recalibration.out.g_9_outputFileTxt02_g_13
g_9_bamFile10_g_13 = Base_Recalibration.out.g_9_bamFile10_g_13


BQSR(g_9_bamFile10_g_13,g_10_genomeDict01_g_13,g_9_outputFileTxt02_g_13)
g_13_bamFile00_g_17 = BQSR.out.g_13_bamFile00_g_17


prepare_comparisons(g_13_bamFile00_g_17.collect(),g_18_1_g_17,g_19_2_g_17)
g_17_inputDir00_g_21 = prepare_comparisons.out.g_17_inputDir00_g_21


mutect2(g_17_inputDir00_g_21.flatten(),g_10_genomeDict01_g_21)
g_21_vcfFile00_g_23 = mutect2.out.g_21_vcfFile00_g_23


vep(g_21_vcfFile00_g_23,g_2_1_g_23,g_24_2_g_23,g_27_3_g_23)
g_23_vcfFile00 = vep.out.g_23_vcfFile00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
