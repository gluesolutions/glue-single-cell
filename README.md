# single cell data
Support for single cell data in glue.

Contents:
1. A data loader to import .h5ad or .loom files into glue as a set of linked glue Data Objects
2. A custom wrapper for the X array in an AnnData file-type to allow a subset of out-of-memory operations 
3. A sample menubar plugin to enable running Scanpy differential gene expression on glue user-defined subsets
4. A sample menubar plugin to enable calculating summary statistics on cells based on a subset of genes

# Installation

If installing into a new virtual environment:

`conda create -n single-cell python=3.9`
`conda activate single-cell`
`conda install -c glueviz glueviz`
`pip install git+https://github.com/gluesolutions/glue-single-cell.git`

Otherwise, the (non-glue) requirements for this package are just `scanpy` and `anndata` so in an environment already running glue the following is all that is required:

`pip install git+https://github.com/gluesolutions/glue-single-cell.git`

# Usage

The Scanpy differential gene expression plugin requires some custom actions to run at startup, which can be achieved by running glue using this command:

`glue -v --startup=setup_anndata`

## Example Doing Scanpy differential gene expression on glue user-defined subsets.
![Differential-Gene-Expression-Overview-Ann](https://user-images.githubusercontent.com/3639698/160698692-258365f1-e9f1-488b-9b92-24b1a0429c47.png)

