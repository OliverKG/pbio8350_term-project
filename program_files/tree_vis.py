import sys
from itolapi import Itol
from itolapi import ItolExport
from pathlib import Path

itol_uploader = Itol()
itol_uploader.add_file(Path(sys.argv[1]))
itol_uploader.params["treeName"] = "tree1"

success = itol_uploader.upload()

if(success):
    itol_exporter = ItolExport()
    itol_exporter.add_export_param_value('tree',itol_uploader.comm.tree_id)
    itol_exporter.add_export_param_value('format',"png")
    itol_exporter.export(Path(sys.argv[2]))
else: print("failure")
