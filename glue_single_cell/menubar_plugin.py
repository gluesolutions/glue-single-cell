from glue.config import menubar_plugin

from .qt.diff_gene_exp import DiffGeneExpDialog


@menubar_plugin("Scanpy Differential Gene Expression")
def diff_gene_exp_plugin(session, data_collection):
    DiffGeneExpDialog.create_subset(data_collection,
                                    default=None, parent=None)
    return
