from pathlib import Path
from eMolFragTEMP.src.input import AcquireFiles

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
    argTypes = ["-i","-o","-p","-m","-c"]
    
    #parse through each argument and send them into the handler
    for i in range(len(arguments)):
        if (arguments[i] in argTypes):
            initializer.setOption(arguments[i], arguments[i + 1])

    return True
            
#
# Reads user input and begins configuration
#
def readConfigurationInput(initializer, arguments):
    arg = arguments
    requiredTypes = ["-i", "-o"]

    #if length is 1, then no arguments wer eprovided
    if (len(arg) <= 1):
        print(f"No arguments we're provided")
        return False

    #if length is 2, then only one argument was provided, meaning the only argument is a file
    if (len(arg) == 2):
        arg = readConfigurationFile(arg[1])

    #after a file can be read, check if we have an accurate amount of arguments
    if (len(arg) <= 4):
      print(f"Too few command arguments") 
      return False

    #if there is no "-i" or "-o", throw an error
    if not(all(x in arg for x in requiredTypes)):
      print(f"Every command must iclude '-i' and '-o'")
      return False

    #otherwise read the command line arguments provided
    readCommandLine(initializer, arg)
    return True