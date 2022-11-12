from glue.config import viewer_tool
from glue.viewers.common.tool import Tool

from glue_heatmap.qt.extract_to_matrix import ExtractToMatrixDialog

PARENT_STRAINS = ['AJ','B6','129S1','NOD','NZO','CAST','PWK','WSB']

@viewer_tool
class QTLSubsetTool(Tool):

	icon = 'glue_spawn'
	tool_id = 'qtl:subset'
	action_text = 'Create a new 2D Matrix a Subset of QTL data'
	tool_tip = 'Create a new 2D Matrix a Subset of QTL data'
	shortcut = 'Ctrl+L'

	def __init__(self, viewer):
		super().__init__(viewer)
		self._data_collection = viewer._data
		self.viewer  = viewer
		
	def activate(self):
		"""
		Fired when the toolbar button is activated
		"""
		self.data = self.viewer.state.layers_data[0]
		default_components = []
		for parent_strain in PARENT_STRAINS:
			default_components.append(self.data.id[parent_strain])
		ExtractToMatrixDialog.create_matrix(
										self._data_collection,
										default_data=self.data, 
										default_components=default_components,
										default_row=self.data.id['gene_id'],
										default_col_names="Parent Strain",
										parent=None)

	def close(self):
		if hasattr(self.viewer, 'window_closed'):
			self.viewer.window_closed.disconnect(self._do_close)
		self.viewer = None
