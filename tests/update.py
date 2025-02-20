

from DataPrepare import UpdatePanel
import os


def run_update():

    password = input("Enter admin password: ")
    
    if len(password) == 3:
        
        up = UpdatePanel()
        
        _cwd = os.getcwd()
        up.path_inside = os.path.join(_cwd, 'data')
        up.path_in_inside = os.path.join(_cwd, 'data', 'in_use')
        up.path_tmp = os.path.join(_cwd, 'data', 'tmp')
                
        
        up.update_from_sources(admin_user = password)
        
    else:
        
        raise ValueError("Incorrect password length!")
    


if __name__ == "__main__":
    run_update()