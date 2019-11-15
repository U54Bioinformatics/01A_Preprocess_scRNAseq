# Preprocess 10x data.
betsy_run.py --num_cores 20 --network_png count12.pdf --receipt count14.txt \
    --input TenXFastqFolder --input_file fastq01 \
    --input SampleGroupFile --input_file samp17.xls \
    --output TenXPreprocessing --output_file count11 \
    --mattr sample_names_are_10x_format=yes \

# Preprocess iCell8 data.
GTF=genomes/GENCODE.GRCh38.v29.annotation.gtf
betsy_run.py --num_cores 30 \
  --network_png proc02.pdf --receipt proc03.txt \
  --input ICell8FastqFolder --input_file fastq21 \
  --input ICell8WellListFile --input_file well81.txt \
  --input ICell8SampleFile --input_file samp81.txt \
  --input GTFGeneModel --input_file $GTF \
  --input ReferenceGenome --input_file genomes/GENCODE.GRCh38.p12 \
  --output ICell8Preprocessing --output_file proc01 \
  --dattr ICell8Preprocessing.aligner=star \
  --dattr ICell8Preprocessing.gene_expression_estimator=featurecounts \
  --mattr num_reads_in_subset=50000 \
  --mattr max_files_in_subset=20 \
  --mattr use_shared_memory=yes



# Do copy number analysis.
betsy_run.py --network_png net.pdf --num_cores 20 \
  --input GenericCNVResults --input_file cnv \
  --input GenericCNVModelSelectionFile --input_file model.txt \
  --input ReferenceGenome --input_file Broad.hg19 \
  --input GTFGeneModel --input_file gtf.txt \
  --output CopyNumberAnalysis —output_file cnv.out \
  —dattr FACETSModelSelectionFile.model_selection=adhoc \
  —mattr facets_gbuild=hg19 \
  —mattr cn_header=tcn.em \
  —mattr total_cn_header=tcn.em \
  —mattr minor_cn_header=lcn.em \
  —mattr cn_header2=CNt \
  —mattr total_cn_header2=CNt \
  —mattr minor_cn_header2=B \
  —mattr discard_chrom_with_prefix=GL,JH,KB,KE,NC_,MT



# Do variant calls from FASTQ files with freebayes.

A=Mills_and_1000G_gold_standard
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
  --input FastqFolder --input_file proc21 \
  --input SampleGroupFile --input_file samp41.txt \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr BamFolderChunked.base_quality_recalibrated=yes \
  --dattr VCFFolder.caller=freebayes \
  --dattr VCFFolder.vartype=snp \
  --dattr SimpleVariantMatrix.caller_suite=single \
  --mattr wgs_or_wes=wgs \
  --mattr realign_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr realign_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr recal_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites3=v/dbsnp_138.b37.vcf.gz 

The input files you need:
- “proc21” is a directory that contains all your Fastq files
- “samp41.txt” is a file that you create that tells which fastq file goes with which sample.  It can be a tab-delimited text file, or an excel file.  I have attached an example file.  You can use the same format, but use your data instead.
- “Homo_sapiens_assembly19.fasta” is the reference genome.
- Mills_and_1000G_gold_standard.indels.b37.vcf.gz
 1000G_phase1.indels.b37.vcf.gz
 dbsnp_138.b37.vcf.gz
 These are files that come with the reference genome.  They tell you where the SNPs and indels are.  It is used for indel realignment.

There’s a line:
--mattr wgs_or_wes=wgs
If you are running with exome sequencing, please change to:
--mattr wgs_or_wes=wes


# Do variant calls from BAM files with freebayes.
A=Mills_and_1000G_gold_standard
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
  --input BamFolder --input_file proc21 \
  --dattr BamFolder.aligner=bwa_mem \
  --input SampleGroupFile --input_file samp41.txt \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr BamFolderChunked.base_quality_recalibrated=yes \
  --dattr VCFFolder.caller=freebayes \
  --dattr VCFFolder.vartype=snp \
  --dattr SimpleVariantMatrix.caller_suite=single \
  --mattr wgs_or_wes=wgs \
  --mattr realign_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr realign_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr recal_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites3=v/dbsnp_138.b37.vcf.gz 



# Do variant calls from FASTQ files with freebayes, no indel
# realignment.
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
 --input FastqFolder --input_file proc21 \
 --input SampleGroupFile --input_file samp41.txt \
 --input ReferenceGenome \
 --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
 --output SimpleVariantMatrix --output_file call01.txt \
 --also_save_highest ManyCallerVCFFolders,call05 \
 --dattr BamFolderChunked.base_quality_recalibrated=no \
 --dattr BamFolderChunked.indel_realigned=no \
 --dattr VCFFolder.caller=freebayes \
 --dattr VCFFolder.vartype=snp \
 --dattr SimpleVariantMatrix.caller_suite=single \
 --mattr wgs_or_wes=wgs


# Do variant calls from BAM files with freebayes, no indel
# realignment.
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
--input BamFolder --input_file proc21 \
--dattr BamFolder.aligner=bwa_mem \
--input SampleGroupFile --input_file samp41.txt \
--input ReferenceGenome \
--input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
--output SimpleVariantMatrix --output_file call01.txt \
--also_save_highest ManyCallerVCFFolders,call05 \
--dattr BamFolderChunked.base_quality_recalibrated=no \
--dattr BamFolderChunked.indel_realigned=no \
--dattr VCFFolder.caller=freebayes \
--dattr VCFFolder.vartype=snp \
--dattr SimpleVariantMatrix.caller_suite=single \
--mattr wgs_or_wes=wgs



# Filter mutations for PyClone.
betsy_run.py --network_png pc02.pdf --num_cores 40 \
   --input SimpleVariantMatrix --input_file mut07.txt \
   --dattr SimpleVariantMatrix.with_coverage=yes \
   --input SequenzaResults --input_file cn01 \
   --input SequenzaModelSelectionFile --input_file mod02.xls \
   --output PyCloneMutationsFolder --output_file pc01 \
   --mattr cn_header=CNt \
   --mattr total_cn_header=CNt \
   --mattr minor_cn_header=B \
   --mattr filter_by_min_total_reads=20 \
   --mattr filter_by_min_alt_reads=5 \
   --mattr filter_by_min_vaf=0.05 \
   --mattr max_copynum_for_pyclone=8 \
   --mattr use_only_consistent_cn=yes

# Run PyClone
GTF=genomes/Broad.hg19.RefSeq.NM_only.no_isoforms.170703.gtf
betsy_run.py --network_png pc02.pdf --num_cores 40 \
  --input SimpleVariantMatrix --input_file mut07.txt \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --input SequenzaResults --input_file cn01 \
  --input SequenzaModelSelectionFile --input_file mod02.xls \
  --input GTFGeneModel --input_file $GTF \
  --output PyCloneAnalysis --output_file pc11 \
  --mattr cn_header=CNt \
  --mattr total_cn_header=CNt \
  --mattr minor_cn_header=B \
  --mattr filter_by_min_total_reads=20 \
  --mattr filter_by_min_alt_reads=5 \
  --mattr filter_by_min_vaf=0.05 \
  --mattr max_copynum_for_pyclone=8 \
  --mattr use_only_consistent_cn=yes


# Differential expression analysis.
betsy_run.py --network_png test02.pdf --num_cores 20 \
  --input SignalFile --input_file exp01.txt \
  --dattr SignalFile.preprocess=counts \
  --input SimpleClassFile --input_file class01.xlsx \
  --dattr SimpleClassFile.preprocess=counts \
  --output DiffExpAnalysis --output_file test01 \
  --dattr DiffExpAnalysis.de_algorithm=deseq2 \
  --dattr DiffExpAnalysis.de_comparison=one_vs_others \
  --dattr DiffExpAnalysis.preprocess=counts \
  --dattr SignalFile.filter_and_threshold_genes=yes \
  --mattr keep_genes_expressed_in_perc_samples=5 \
  --mattr p_cutoff=0.05 \
  --mattr plot_max_genes_per_group=10


# Run Sequenza.
betsy_run.py --network_png cn12.pdf --receipt cnv13.txt --num_cores 20 \
   --input BamFolder --input_file bam11 \
   --dattr BamFolder.has_read_groups=yes \
   --dattr BamFolder.sorted=coordinate \
   --input ReferenceGenome --input_file genomes/Broad.hg19 \
   --input NormalCancerFile --input_file nc01.xls \
   --output SequenzaResults --output_file cn11 \
   --mattr sequenza_assembly=hg19 \
   --mattr discard_chrom_with_prefix=GL,JH,KB,KE,NC_,MT



# UMAP on ssGSEA scores.
betsy_run.py --network_png clust12.pdf --num_cores 40 \
   --input SignalFile --input_file gsea22.txt \
   --dattr SignalFile.has_entrez_gene_info=yes \
   --dattr SignalFile.has_human_gene_info=yes \
   --dattr SignalFile.logged=not_applicable \
   --output UMAPAnalysis --output_file clust11 \
   --mattr umap_alternate_pca_genes=10,25 \
   --mattr umap_alternate_pca_dimensions=5,10 \
   --mattr umap_alternate_snn_k=10,25,50 \
   --mattr umap_alternate_snn_resolution=0.6,0.8,1.0 \
   --mattr umap_alternate_umap_dimensions=5,10 \
   --mattr umap_pca_genes=25 \
   --mattr umap_pca_dimensions=5 \
   --mattr umap_snn_k=10 \
   --mattr umap_snn_resolution=0.8 \
   --mattr umap_umap_dimensions=10



# Somatic variant calls from BAM files.
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz
betsy_run.py --num_cores 20 --network_png call02.pdf --receipt call03.txt \
  --input BamFolder --input_file bam01 \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --input SampleGroupFile --input_file samp02.xlsx \
  --input NormalCancerFile --input_file nc01.xlsx \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr SimpleVariantMatrix.duplicates_marked=yes \
  --dattr SimpleVariantMatrix.caller_suite=lance3 \
  --mattr wgs_or_wes=wgs \
  --mattr mutect_dbsnp_vcf=MuTect/dbsnp_132_b37.leftAligned.vcf \
  --mattr mutect_cosmic_vcf=MuTect/b37_cosmic_v54_120711.vcf \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --mattr filter_by_min_total_reads=20 \
  --dattr SimpleVariantMatrix.annotated_with_annovar=yes \
  --mattr annovar_buildver=hg19 \
  --dattr SimpleVariantMatrix.annotated_with_snpeff=yes \
  --mattr snpeff_genome=GRCh37.75 \
  --dattr SimpleVariantMatrix.variant_annotations=cancer_cosmic \
  --mattr cancer_genes_file="008 Cancer Genes/cancer_genes.txt" \
  --mattr cosmic_variants_file="008 Cancer Genes/${COSMIC}" \
  --dattr SimpleVariantMatrix.with_coverage=yes 


