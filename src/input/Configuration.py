from pathlib import Path
from eMolFragTEMP.src.input import AcquireFiles
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.utilities import logging

def checkRequirements(arg):
  if ((arg.i == None) or (arg.o == None)):
    logging.logger.error(f"Every command must iclude '-i' and '-o'")
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
        logging.logger.error(f"File {path.name} is empty")
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
      logging.logger.error(f"Arguments were empty")
      return None

    initializer.setOptions(arguments)

    return initializer
            
#
# Reads user input and begins configuration
#
def readConfigurationInput(initializer, arguments):
    arg = arguments
    
    if (arg.i == None):
      logging.logger.error(f"No input was provided")
      return None
    inPath = Path(arg.i)
    try:
        if (inPath.suffix == ".txt"):
          retArg = readConfigurationFile(arg.i)
          arg = utilities.createParser(retArg)
    except:
        pass
    
    if not(checkRequirements(arg)):
      return None

    #otherwise read the command line arguments provided
    readCommandLine(initializer, arg)
    return initializer
