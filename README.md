# Somatic Variant Calling Pipeline

#### Identify somatic short variants (SNVs and Indels) in a matched tumor and normal sample. (<a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels" target="_blank" rel="noopener noreferrer"><span style="color:blue">documentation</span></a>)

##### **Module Input:**

1) FASTQ files for each sample (both tumor and matching normal).

2)  A "groups" file which contains sample information. **Requirements:** 
  * Tab or comma separated (TSV, CSV) with a header.
  * **sample_name** and **group** columns are required.
  * Each sample_name must appear as a column in the header of the counts file.
  * Additional columns can be included which contain metadata about each sample.
  * Values in **sample_name** column cannot contain spaces.

3) A "comparison" file which specifies which groups to compare in SNP/indel analysis. **Requirements**:  
  * Tab or comma separated (TSV, CSV) with a header.
  * column names are exactly:
  
   | controls | treats | names | 
	 
  * The values in treats and controls columns must exist in the **group** column from the groups file.
  * The names column will control the prefix of output files.
  * The **names** column cannot contain characters that are forbidden in file names (spaces, %)

4) A list of known SNPs and known Indels in vcf.gz form with accompanying indexes.

**System inputs for custom genome build:**
* GTF File for Variant Effect Predictor (VEP): https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff



