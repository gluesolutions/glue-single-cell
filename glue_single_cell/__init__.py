def setup():
    from .anndata_factory import read_anndata # noqa
    from .anndata_factory import setup_anndata # noqa
    from .menubar_plugin import diff_gene_exp_plugin  # noqa
    from .menubar_plugin import pca_subset_exp_plugin  # noqa
    from .menubar_plugin import enrichrpy_plugin # noqa
    from .qtl_viewer.viewer import QTLViewer # noqa

    from glue.config import qt_client
    qt_client.add(QTLViewer)

    # Add colormaps we can use when we use the pca_subset_exp_plugin 
    # to maps gene expression to a summary over cells. A summary generated 
    # from a red subset can then be displayed using the 'Reds' colormap
    # This does not happen automatically... at least not yet

    from glue.config import colormaps
    import matplotlib.cm as cm
    colormaps.add('Reds', cm.Reds)
    colormaps.add('Greens', cm.Greens)
    colormaps.add('Blues', cm.Blues)

    # Add some categorical colormaps

    colormaps.add('Pastel1', cm.Pastel1)
    colormaps.add('Paied', cm.Paired)
    colormaps.add('tab20', cm.tab20)