# singlecell-scripts
My single cell RNA and RNA/ATAC scripts for people who are new to single cell analysis

I would recommend running through scmultiome.R if you're new because it has the most detailed explanations and only deals with one sample. You can ask me for a sample dataset. Below is my recommended workflow, but if you only have one sample and don't need to integrate you can skip objectintegration.R.

<img width="403" height="131" alt="image" src="https://github.com/user-attachments/assets/545ff357-76b9-4a01-a708-8e30024de8c2" />


## importparse.R
This is for if you're working with Parse sequencing technology and have downloaded the matrices after the initial alignment to analyze yourself. 

## objectintegration.R
If you have more than one sample that you need to integrate, use objectintegration.R. This includes scripts for integrating one or more scRNA samples as well as one or more scMultiome (RNA+ATAC) samples. Also importantly, it gives you information on how to subset a portion of your integrated object and renormalize. This method also works for subsetting a non-integrated object. 

## pseudotimetrajectory.R
Contains scripts for Monocle3 and Slingshot.

## misc.R
Contains scripts for:
* reading in different file types to make a Seurat object (.mtx, .hdf5, .h5ad)
* transferring metadata identities between two objects based on cell barcodes (both objects must have some of the same cells)
* adding metadata categories based on expression of a single gene
* visualizing PCA dimensions
* specifying colors for any given plot
* defining ordered levels of cluster identities
* making a variability plot using HVGs (high variariability genes)
* visualizing module scoring based on a gene list
* changing ENSEMBL IDs to gene symbols.
