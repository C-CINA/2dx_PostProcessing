from argparse import ArgumentParser

from SystemUtils import module_path


class Input:
    '''
    An input class to hold the data  for input.
    The class contains following elements:\n
    name   : contains variable name for positioned (compulsory) variables
                  and variable name STARTTING WITH - for optional arguments
    identifier : one letter identifier only if optional
    default    : Default value if optional
    type       : type of the variable (int/float)
    help       : help text for the variable
    '''
    
    def __init__(self, name, identifier=None, default=None, inptype=None, helpx=""):
        self.name = name
        self.identifier = identifier
        self.default = default
        self.type = inptype
        self.helpx = helpx

    
class InputParser:
    '''
    A class used to parse the inputs. 
    Possible inputs are read from a file.
    '''
    inputs = []

    def __init__(self, input_list_file=module_path() + "/inputs.default_value", args=None):
        self._add_input(input_list_file)
        self._init_parser()
        self.arguments = self._parse(args)
        
    def _add_input(self, input_list_file):
        '''
        Reads in the possible inputs from a file. The file should have following format:
        COLUMNS:
        -----------------------------------------------------------
        VARIABLE_NAME  IDENTIFIER  DEFAULT_VALUE  TYPE  HELP
        -----------------------------------------------------------
        *use None for default and identifier if input is compulsory
        *Put None whenever a field is not required
        *Supported types: int, str, float, None
        '''
        f = open(input_list_file)
        for line in f.readlines():
            if not line.startswith("#") and line.strip():
                
                try:
                    [var, idx, default_r, type_r, help_r] = ' '.join(line.split()).split(' ', 4)
                except:
                    print 'WARNING: Error in reading file specifying inputs, Please check:'
                    print '{}\n' .format(line)
                    continue
                
                if help_r is not 'None':
                    helpx = help_r
                else:
                    helpx = ""
                  
                if default_r == 'None':
                    name = var
                    identifier = None
                    default = None
                else:
                    name = var
                    identifier = idx
                    default = default_r
                    
                if type_r == 'int':
                    inptype = int
                    default = int(default)
                elif type_r == 'str':
                    inptype = str
                elif type_r == 'float':
                    inptype = float
                    default = float(default)
                else:
                    inptype = None
            
                self.inputs.append(Input(name, identifier, default, inptype, helpx))
        
        f.close()          

    def _init_parser(self):
        description = 'A program which takes a hkl file as input ' \
                      'does post processing and outputs final the 3D-Image(.mrc) ' 
        
        self.parser = ArgumentParser(description=description)
        
        for inputx in self.inputs:
            if inputx.identifier == None:
                flag = inputx.name
                self.parser.add_argument(flag, help=inputx.helpx)
            else:
                self.parser.add_argument('-' + inputx.identifier, '--' + inputx.name,
                                         default=inputx.default,
                                         type=inputx.type,
                                         help=inputx.helpx)

    def _parse(self, args):
        return self.parser.parse_args(args)
    
    def print_inputs(self):
        '''
        Prints all the inputs used by the system
        '''
        print 'Using the following parameters:'
        for attr in self.inputs:
            print '{:20s}: {}' .format(attr.name, getattr(self.arguments, attr.name))      
        print '\n'
            
            
