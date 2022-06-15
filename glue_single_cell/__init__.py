def setup():
    from .anndata_factory import read_anndata # noqa
    from .anndata_factory import setup_anndata # noqa
    from .menubar_plugin import diff_gene_exp_plugin  # noqa
    from .menubar_plugin import pca_subset_exp_plugin  # noqa
    from .qtl_viewer.viewer import QTLViewer # noqa
    
    from glue.config import qt_client
    qt_client.add(QTLViewer)