
from Options import Options
from pathlib import Path

def acquireFiles(initializer):

    folderPath = Path(initializer.INPUT_PATH)

    if not folderPath.exists():
        print(f'Input path {initializer.INPUT_PATH} does not exist.')
        return
    
    # initialize two lists; one with the full path of a file and one with the file names
    files = []
    filenames = []
    bad_files = []
    for file in folderPath.iterdir():
    
        # file extension will help filter bad data
        extension = file.suffix

        if extension not in constants.ACCEPTED_FORMATS:
            bad_files.append(file)

        else:
            filenames.append(file.name)
            files.append(Path(initializer.INPUT_PATH + "/" + file.name))
            
    if bad_files:
        print(f'[Error] emolFrag 2.0 only accepts the following formats {', '.join(constants.ACCEPTED_FORMATS)}')
        print(f'The following files will be ignored:\n\t{"\n\t".join(bad_files)}')
            
    return files