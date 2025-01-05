# manualy data prepare and adjustment for validated version of GEDS data base

from DataPrepare import UpdatePanel

up = UpdatePanel()

up.check_last_update()


up.update_from_sources(admin_user = None)


