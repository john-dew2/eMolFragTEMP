
from eMolFragTEMP.src.input import Options
from pathlib import Path
from eMolFragTEMP.src.utilities import constants

#
# Given the input path conatining molecules, make a list of all the file paths
#
def acquireMoleculeFiles(initializer):

    folderPath = Path(initializer.INPUT_PATH)

    if not folderPath.exists():
        print(f'Input path {initializer.INPUT_PATH} does not exist.')
        return None
    
    # initialize two lists; one with the full path of a file and one with the file names
    files = []
    filenames = []
    bad_files = []
    for current_file in folderPath.iterdir():
    
        # file extension will help filter bad data
        extension = current_file.suffix

        if extension not in constants.ACCEPTED_FORMATS:
            bad_files.append(current_file)

        else:
            filenames.append(current_file.name)
            files.append(Path(initializer.INPUT_PATH + "/" + current_file.name))
            
    if bad_files:
        print(f"[Error] emolFrag 2.0 only accepts the following formats {', '.join(constants.ACCEPTED_FORMATS)}")
        print(f"The following files will be ignored:\n\t{'  '.join(bad_files)}")
            
    return files

#
# Given a configuration file, return the file path
#
def acquireConfigurationFile(file):
    filePath = Path(file)

    if not filePath.exists():
        print(f'Input path {file} does not exist.')
        return None
    
    return filePath