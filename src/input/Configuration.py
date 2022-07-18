from pathlib import Path
from eMolFragTEMP.src.input import AcquireFiles
from eMolFragTEMP.unittests import utilities

def createParser():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='eMolFrag2')
    parser.add_argument("-i",
    type=str,
    help='Set the input path')
  
    parser.add_argument("-o",
    type=str,
    help='Set the output path')
  
    parser.add_argument("-m",
    type=int, choices=range(0,3),
    help='Set the execution type')
  
    parser.add_argument("-c",
    type=int, choices=range(0,3),
    help='Set the output type')

    return parser

def checkRequirements(arg):
  if ((arg.i == None) or (arg.o == None)):
    utilities.emit(0, f"Every command must iclude '-i' and '-o'")
    return False
  return True
#
# Acquires and reads a configuration file and returns the command line arguments nested in the file
#
def readConfigurationFile(config_file):
    #acquire path to folder
    path = AcquireFiles.acquireConfigurationFile(config_file)
    retString = ""
    
    #read the lines
    with open(path) as f:
        lines = f.readlines()

    if (len(lines) <= 0):
        utilities.emit(0, f"file {path.name} is empty")
        return []
    
    #concatenate the contents and ignore comments
    position = 0
    for line in lines:
        position = line.find("#")
        if (position > -1):
            retString += line[:position]
        else:
            retString += line
    
    #fragmentize the string into a list of command arguments
    return retString.split(" ")

#
# Reads a command line and set options to correct values
#
def readCommandLine(initializer, arguments):
    if (arguments == None):
      utilities.emit(0, f"Arguments were empty")
      return None

    initializer.setOptions(arguments)

    return initializer
            
#
# Reads user input and begins configuration
#
def readConfigurationInput(initializer, arguments):
    arg = arguments

    if (arg.i == None):
      utilities.emit(0, f"No input was provided")
      return None

    inPath = Path(arg.i)
    try:
        if (inPath.suffix == ".txt"):
          retArg = readConfigurationFile(arg.i)
          parser = createParser()
          arg = parser.parse_args(retArg)
    except:
        pass
    
    if not(checkRequirements(arg)):
      return None

    #otherwise read the command line arguments provided
    readCommandLine(initializer, arg)
    return initializer
