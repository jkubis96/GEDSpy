
from data_prepare import Donwload

dw = Donwload()

#
ref = dw.download_ref()

del ref

#
rns_seq = dw.download_rns_seq()

del rns_seq

#
IntAct = dw.download_IntAct_data()

del IntAct

"""
This method retrieves data from the IntAct database 

Returns:
   dict: IntAct database file in dictionary format
   
"""
     
   
#
Viral = dw.download_viral_deiseases()

del Viral

"""
This method retrieves data from the ViMIC database 

Returns:
   dict: ViMIC database file in dictionary format
   
"""


#
HPA = dw.download_HPA()

del HPA

"""
This method retrieves data from the Human Protein Atlas database 
   
Returns:
   dict: Human Protein Atlas database file in dictionary format
   
"""


#
string = dw.download_string()

del string

"""
This method retrieves data from the STRING database 
   
Returns:
   dict: STRING database file in dictionary format
   
"""


#
kegg = dw.download_kegg()

del kegg

"""
This method retrieves data from the KEGG database 

Returns:
   dict: KEGG database file in dictionary format
   
"""


#
reactome = dw.download_reactome()

del reactome

"""
This method retrieves data from the REACTOME database 
   
Returns:
   dict: REACTOME database file in dictionary format
   
"""


#
go = dw.download_go_term()

del go

"""
This method retrieves data from the GO-TERM database 

Returns:
   dict: GO-TERM database file in dictionary format
   
"""


#
cp = dw.download_cell_phone()

del cp

"""
This method retrieves data from the CellPhone database 
   
Returns:
   dict: CellPhone database file in dictionary format
   
"""


#
ct = dw.download_cell_talk()

del ct

"""
This method retrieves data from the CellTalk database 
   
Returns:
   dict: CellTalk database file in dictionary format
   
"""



dw.update_downloading(password = 'JBS')
 


dw.check_last_update()


"""
This method checks the last udate of data used in this library
   
Returns:
   date: Date of last update
   
"""



from data_prepare import DataAdjustment

da = DataAdjustment()


da.update_to_data()


# zip data

da.ZIP()


# move zip to different directory
import os

da.get_ZIP(path = os.getcwd())



