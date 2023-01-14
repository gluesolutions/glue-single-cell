import os
import subprocess
import re
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt
from glue.utils.qt import load_ui
from qtpy.QtWidgets import QDialog, QVBoxLayout
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5 import QtCore
from glue.core.state_objects import State
from echo import SelectionCallbackProperty, CallbackProperty
from glue.core.data_combo_helper import DataCollectionComboHelper, ComponentIDComboHelper, ComboHelper, ManualDataComboHelper

__all__ = ['GenomeTrackViewerDialog', 'GenomeTrackViewerState']


class GenomeTrackViewerState(State):
    data = SelectionCallbackProperty() 
    subset = SelectionCallbackProperty()

    def __init__(self, data_collection):
        
        super().__init__()

        self.data_collection = data_collection
        self.data_helper = DataCollectionComboHelper(self, 'data', data_collection)
        self.subset_helper = ComboHelper(self, 'subset')

        def display_func_label(subset_group):
            return subset_group.label
        
        self.add_callback('data', self._on_data_change)
        self._on_data_change()
        
        self.subset_helper.choices = data_collection.subset_groups
        try:
            self.subset_helper.selection = data_collection.subset_groups[0]
        except IndexError:
            pass
        self.subset_helper.display = display_func_label
        
    def _on_data_change(self, *args, **kwargs):
        pass


class GenomeTrackViewerDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super().__init__(parent=parent)

        self.state = GenomeTrackViewerState(collect)

        self.ui = load_ui('ucsc_genome_track_viewer.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self, do_dialog=True):
        """
        Extract the relevant piece of the datafile to a BED file
        """
        if self.state.subset is not None:
            for subset in self.state.subset.subsets:
                if subset.data == self.state.data:
                    gene_subset = subset
        gene_names = gene_subset[self.state.data.id['gene_id']]
        gene_starts = gene_subset[self.state.data.id['gene_start']]
        gene_ends = gene_subset[self.state.data.id['gene_end']]
        browser_lines = "browser position chr5:137000000-150000000"
        browser_lines_2 = "browser hide all"
        browser_lines_3 = "browser dense knownGene"
        track_lines = "track name=coords description='Gene Subset from glue' color=96,39,254"
        header_lines = "#chrom chromStart chromEnd name"
        f = open('temp.bed', 'at', encoding='utf-8')
        for header_string in [browser_lines, browser_lines_2, browser_lines_3, track_lines, header_lines]:
            f.write(header_string)
            f.write(os.linesep)
        for name, start, end in zip(gene_names, gene_starts, gene_ends):
            if start < 151758149 and end < 151758149:
                f.write(f"chr5\t{start}\t{end}\t{name}")
                f.write(os.linesep)
        f.close()
        
        p1 = subprocess.run(["curl", "-s", "-F", "db=mm39", "-F", 'hgct_customText=@temp.bed', "http://genome.ucsc.edu/cgi-bin/hgCustom"],capture_output=True)        
        for line in p1.stdout.decode('utf-8').split('\n'):
            if 'hgsid=' in line:
                hgsid = re.search('hgsid=(.*)"', line).group(1)
                break
                print(hgsid)
        start_url = f'http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid={hgsid}&position=chr5:137000000-150000000'
        self.web = QWebEngineView()
        self.web.setUrl(QtCore.QUrl(start_url))
        self.web.resize(1000,800)
        self.web.show()
        os.remove('temp.bed')
        #self.layout = QVBoxLayout()
        #self.layout.addWidget(self.web)
        #self.layout.setContentsMargins(3, 3, 3, 3)
        
        #self.setLayout(self.layout)

    @classmethod
    def display(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
