# single cell data
Support for analysis of single cell data in glue.

Contents:
1. A data loader/dialog to import .h5ad or .loom files into glue as a set of linked glue Data objects
2. A custom wrapper for the X array in an AnnData file-type to allow a subset of out-of-memory operations 
3. A menubar plugin to enable running Scanpy differential gene expression on glue user-defined subsets
4. A menubar plugin to enable calculating summary statistics on cells based on a subset of genes
5. A menubar plugin to add KEGG pathways to a gene list
6. A custom viewer for QTL data, including chromosome-labelling and filtering by LOD
7. A viewer for small multiples/facets for examination of different experimental factors.

# Installation

## Recommended

We recommend installing this plug-in into a new virtual environment:

`conda create -n single-cell`

`conda activate single-cell`

`conda install -c glueviz glueviz`

`pip install git+https://github.com/gluesolutions/glue-single-cell.git`

## Manual

If installing into an exsiting glue environment the additional requirements for this package are `scanpy`,
`anndata`, and `enrichrpy` which can generally be installed via simply as:

`pip install git+https://github.com/gluesolutions/glue-single-cell.git`

# Usage

The Scanpy differential gene expression plugin requires some custom actions to run at startup, which can be achieved by running glue using this command:

`glue -v --startup=setup_anndata`


## Data Structure

This plug-in maps the [AnnData schema](https://anndata.readthedocs.io/en/latest/) for storing single-cell data into linked [glue Data objects](http://docs.glueviz.org/en/stable/python_guide/data_tutorial.html) as follows:

![DataStructureSchematic](https://user-images.githubusercontent.com/3639698/164315869-935163b1-2503-4e12-8166-3978da8dbe0c.png)

This structure allows for using all the normal glue viewers and functionality on the tabular obs and var arrays while retaining access to the X array for calculations that require it.

## Example Doing Scanpy differential gene expression on glue user-defined subsets.
![Differential-Gene-Expression-Overview-Ann](https://user-images.githubusercontent.com/3639698/160698692-258365f1-e9f1-488b-9b92-24b1a0429c47.png)

