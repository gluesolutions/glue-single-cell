from glue.config import menubar_plugin

from .qt.diff_gene_exp import DiffGeneExpDialog
from .qt.pca_subset import PCASubsetDialog
from .qt.gsea import GSEApyDialog


@menubar_plugin("Scanpy Differential Gene Expression")
def diff_gene_exp_plugin(session, data_collection):
    DiffGeneExpDialog.create_subset(data_collection,
                                    default=None, parent=None)
    return

@menubar_plugin("Calculate Summary Over Gene Subset")
def pca_subset_exp_plugin(session, data_collection):
    PCASubsetDialog.summarize(data_collection,
                                    default=None, parent=None)
    return

@menubar_plugin("Enrich Gene Set via Enrichr")
def gseapy_plugin(session, data_collection):
    GSEApyDialog.enrich(data_collection,
                            default=None, parent=None)
    return
