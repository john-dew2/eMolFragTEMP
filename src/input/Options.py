#Options.py
from eMolFragTEMP.src.utilities import logging


#input_path
#output_path
#parallel_cores_used
#outut_type:
    # if 0 then full_process
    # if 1 then chop_only
    # if 2 then chop_and_remove
#output_format:
    #if 0 then traditional_format
    #if 1 then different_files and remove_log_files
    #if 2 then different_files and !remove_log_files
    
class Options:
    def __init__(self):
        self.INPUT_PATH = ""
        self.OUTPUT_PATH = ""
        
        self.INDIVIDUAL = False
        self.ALL = True
        
        self.UNIQUE = False
        self.REDUNDANT = True
        

    #
    # Prints the current preferences for each option available
    #
    def printOptions(self):
        logging.logger.info(f"\nInput Path:  {self.INPUT_PATH}\nOutput Path: {self.OUTPUT_PATH}\nOutput Format:\n    Individual Files:  {self.INDIVIDUAL}\n    All Files:         {self.ALL}\nExecution Type:\n    Remove Redundacy:  {self.UNIQUE}\n    No Removal:        {self.REDUNDANT}")
    #
    # Given an options preference, function sets an option to its correct value
    #
    def setOption(self, argType, option):
    
        if (argType == "-i"):   
            self.INPUT_PATH = option
            
        elif (argType == "-o"):
            self.OUTPUT_PATH = option
            
        elif (argType == "-u"):
            self.INDIVIDUAL = option
            self.ALL = not(option)
        
        elif (argType == "-indiv"):
            self.UNIQUE = option
            self.REDUNDANT = not(option)
        else:
            logging.logger.error(f"{argType} does not exist")

    def setOptions(self, args):
      try:
        self.setOption("-i", args.i)
      except:
        pass
      try:
        self.setOption("-o", args.o)
      except:
        pass
      try:
        self.setOption("-u", args.u)
      except:
        pass
      try:
        self.setOption("-indiv", args.indiv)
      except:
       pass
      #self.printOptions()