import os
import sys
from time import gmtime, strftime

def file_present(file):
    '''
    Checks if a file is present or not with the input path
    '''
    return os.path.isfile(file)

def create_dir(dir_path):
    '''
    Creates a directory if not exist
    '''
    dir_path = dir_path + "/" 
    d = os.path.dirname(dir_path)
    if not os.path.exists(d):
        os.makedirs(d) 

def get_file_name(path):
    '''
    Returns the filename or end directory name of the <path>
    '''
    return os.path.basename(path)
 
def get_dir_path(file):
    '''
    Returns the directory path of the given file
    '''
    return os.path.dirname(os.path.realpath(file))
 
def get_current_time():
    '''
    Returns the current time in the format: YYYY-mm-dd_HH-MM-SS
    '''
    return strftime("%Y-%m-%d_%H-%M-%S", gmtime())

def we_are_frozen():
    '''
    All of the modules are built-in to the interpreter, e.g., by py2exe
    '''
    return hasattr(sys, "frozen")

def module_path():
    '''
    Returns the physical path of the module in consideration
    '''
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))