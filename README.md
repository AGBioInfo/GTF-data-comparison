# GTF-data-comparison

(see Project runbook for result report)

Compare, is ready to run perl script. It takes 2 GTF files as an input and give2 PDF files as an output. One PDF is ‘Plots.pdf’ which have venn diagrams, histograms and some bar plots. Other PDF is ‘DetailedPlots.pdf’ which have more detailed multiple barplots. For data extraction and filtering perl hashes and awk utility is used. To find overlaps intersectBed, one of the tools from BEDTools is used, which is set to find overlap of minimum 50% of base pairs of the shortest feature of 2 files. Here feature is referred to genome features(exon, transcript, gene, CDS,UTR etc). A separate R script is written to generate different type of plots and print them all in one PDF file. This R script is called inside compare perl script, thus no need to run Rscript separately.
