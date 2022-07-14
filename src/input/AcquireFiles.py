
from eMolFragTEMP.src.input import Options
from pathlib import Path
from eMolFragTEMP.src.utilities import constants

#
# Given the input path conatining molecules, make a list of all the file paths
#
def acquireMoleculeFiles(initializer):

  folderPath = Path(initializer.INPUT_PATH)
  
  #if the folder path doesnt exist, exit processing
  if not folderPath.exists():
      print(f'Input path {initializer.INPUT_PATH} does not exist.')
      return None
 
  files = []
  bad_files = []
  #grab every molecule in a folder
  for current_file in folderPath.iterdir():
  
    #if the file extension is not a supportedd format, add the file to the bad file list, otherwise add it to the file list
    extension = current_file.suffix
    if extension not in constants.ACCEPTED_FORMATS:
        bad_files.append(current_file)
    else:
        files.append(Path(initializer.INPUT_PATH + "/" + current_file.name))
  
  #if there were any bad files, print an error
  if (len(bad_files) > 0):
      print(f"[Error] emolFrag 2.0 only accepts the following formats {', '.join(constants.ACCEPTED_FORMATS)}")
      #print(f"The following files will be ignored: {', '.join(bad_files)}")
  
  return files

#
# Given a configuration file, return the file path
#
def acquireConfigurationFile(usr_file):
    filePath = Path(usr_file)

    #if the folder path doesnt exist, exit processing
    if not filePath.exists():
        print(f'Input path {usr_file} does not exist.')
        return None
    
    return filePath