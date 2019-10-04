  betsy_run.py --num_cores 20 --network_png count12.pdf --receipt count14.txt \
    --input TenXFastqFolder --input_file fastq01 \
    --input SampleGroupFile --input_file samp17.xls \
    --output TenXPreprocessing --output_file count11 \
    --mattr sample_names_are_10x_format=yes \



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
