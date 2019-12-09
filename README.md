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


# Preprocess for Fluidigm signal data.
## Make a fluidigm sample file.
 betsy_run.py \
   --input FluidigmFastqFolder --input_file fastq01 \
   --output UnvalidatedFluidigmSampleFile --output_file samp01.xls


 GTF=genomes/Broad.hg19.RefSeq.NM_only.no_isoforms.170703.gtf
 betsy_run.py --num_cores 40 --network_png fastq52.pdf \
   --input FluidigmFastqFolder --input_file fastq01 \
   --input FluidigmSampleFile --input_file samp01.xls \
   --input GTFGeneModel --input_file $GTF \
   --input ReferenceGenome --input_file genomes/Broad.hg19 \
   --output FeatureCountsSummary --output_file fastq51.xls \
   --dattr FeatureCountsSummary.aligner=star \
   --mattr use_shared_memory=yes


# Preprocess for Fluidigm signal reads (already demultiplexed).
 GTF=genomes/gencode.vM14.annotation.gtf
 betsy_run.py --num_cores 20 --network_png proc12.pdf --receipt proc13.txt \
   --input FastqFolder --input_file fastq21 \
   --input SampleGroupFile --input_file samp02.xls \
   --input GTFGeneModel --input_file $GTF \
   --input ReferenceGenome --input_file genomes/GENCODE.GRCm38 \
   --output RNASeqSignalFile --output_file exp01.txt \
   --dattr RNASeqSignalFile.preprocess=counts \
   --dattr RNASeqSignalFile.aligner=star \
   --dattr RNASeqSignalFile.gene_expression_estimator=featurecounts
   
## QC on the data with:
 betsy_run.py --num_cores 20 --network_png meta02.pdf \
   --input SignalFile --input_file count01.txt \
   --dattr SignalFile.preprocess=counts \
   --output scRNAMetadata --output_file meta01



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


