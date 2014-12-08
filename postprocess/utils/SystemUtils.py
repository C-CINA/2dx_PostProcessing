import os
import sys
from time import gmtime, strftime

# Creates a directory if not exist
def create_dir(dir_path):
    dir_path = dir_path + "/" 
    d = os.path.dirname(dir_path)
    if not os.path.exists(d):
        os.makedirs(d)
    return    

# Returns the filename or end directory name of the <path>
def get_file_name(path):
    return os.path.basename(path)

# Returns the directory path of the given file
def get_dir_path(file):
    return os.path.dirname(os.path.realpath(file))

# Returns the current time
def get_current_time():
    return strftime("%Y-%m-%d_%H-%M-%S", gmtime())

def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def module_path():
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))