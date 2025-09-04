R255X E488QD Data Analysis Project
Anthony D. Surkov
Last updated: 04-28-2025


Repository organization:

data/
	archive/					Holds previously-generated and deprecated datasets.
	R255X_E488QD_casey_output.csv			Pre-processed data fed through the Casey pipeline.
	R255X_E488QD_clean.csv				Casey pipeline data fed through a cleaning step to isolate N10, GAA, GGA, edit%, n, n_other.
	R255X_E488QD_modeling.csv			Holds one-hot encoded version of N10 sequences, along with everything in R255X_E488QD_clean.csv.
	R255X_E488QD_preprocesed_SE.fastq(.gz)  	Pre-processed data.
eda/
	reads_diagnostics.R				Supports modeling number of reads per sequence as a negative binomial distribution.
	reads_visualization.R				Visualizes distribution of reads.
	basic_statistics.R				Catch-all program to play with basic EMERGe statistics: read depth, editing frequencies, NA presence, etc.
	pick_n0.R					Exploration into minimum number of reads recommended based on NB distribution found in diagnostics_reads.R.
	edits_visualization.R				Visualizes distribution of editing in UMAP space.

img/
	reads_visualization/				Holds plots associated with diagnostic_reads.R.
	exploratory/					Holds plots associated with exploratory data analysis.
	edits_visualization/					Holds plots associated with the animation of differently-filtered UMAP plots.
integration/
	casey_pipeline.R				Generalized form of the original Casey pipeline. Converts HTStream pre-processed data (R255X_E488QD_preprocesed...) into Casey-output data
	clean_data.R					Converts Casey pipeline data into cleaned data.
	one_hot_encoding.R				Converts cleaned data into one-hot-encoded sequence data.
modeling/
	linear_pairwise_glms.R				Models EMERGe by looking at all single-nucleotide and pairwise-nucleotides as features.
	random_forest.R					Identifies most important features using a random forest.
nd/
	data/	
			archive/			Holds old R255X_E488QD datasets.
			R255X_E488QD_fulldata.csv	Full processed data (n >= 10, edit > 0)
			R255X_E488QD_top_100.csv	Top 100 editors of R255X_E488QD EMERGe screen.
	nd_data_gen.R					Generates R255X_E488QD_fulldata.csv and R255X_E488QD_top_100.csv.
	nd_diagnostics.R				Attempts to identify cause of difference between old AC and Top 35 Editor data and recreated dataset(s).
	nd_validation.R					Compares newly-generated dataset to old AC and top-35 data.