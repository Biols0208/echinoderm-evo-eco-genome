## PGLS Analysis
Rscript PGLS.R analysis_config.yaml

## phylopath Analysis
Rscript phylopath.R

## plot 

### Creates scatter plots with linear regression lines and 95% CI
Rscript plot_lm.R -h
#Rscript plot_lm.R --input genome_Phenotype.csv --xcol DNA --ycol genomesize --pgls true --tree-file ref.tree --species-col Species --outfile output.pdf --pgls-correlation auto --width 13 --height 6.5

###  Creates box plots or violin plots with overlaid points colored by groups
Rscript plot_box.R -h
#Rscript plot_box.R --input genome_Phenotype.csv --xcol Reproductive_Mode --ycol genomesize --groupcol group --outfile test.pdf --theme-style bw --density-plot --effect-size --export-stats
