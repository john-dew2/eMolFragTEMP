
import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.input import AcquireFiles, Configuration, Options
from eMolFragTEMP.src.utilities import logging

usr_dir = Path.cwd()
dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")
config = usr_dir.joinpath("/content/eMolFragTEMP/unittests/data/config-files")

def runAcquireMoleculeFilesTests():
    
    testPaths = [mol2, smi]
    #Tests if mol2 and smi are taken
    for filePath in testPaths:
        runAcquireMoleculeFiles(f"-i {filePath} -o output/".split(" "), 5)
    #Will recognize the files as bad
    runAcquireMoleculeFiles(f"-i {sdf} -o output/".split(" "), 0)
    #File does not exist
    runAcquireMoleculeFiles(f"-i directory/doesnt/exist -o output/".split(" "), 0)
        
def runAcquireMoleculeFiles(arguments, expec):
    arg = utilities.createParser(arguments) 
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arg)
    files = AcquireFiles.acquireMoleculeFiles(initializer)
    if (files == None): files = []
    assert len(files) == expec
    
def runAcquireConfigurationFileTests():
    runAcquireConfigurationFile(config.joinpath("comment1.txt"), config.joinpath("comment1.txt"))
    runAcquireConfigurationFile(config.joinpath("normal1.txt"), config.joinpath("normal1.txt"))
    runAcquireConfigurationFile(mol2.joinpath("no.txt"), None)

def runAcquireConfigurationFile(org_file, expec):
    processed_file = AcquireFiles.acquireConfigurationFile(org_file)
    assert processed_file == expec

#
#
# Unit test initiation functionality
#
#
def runtest1(func):

    try:
       func()
       return True

    except:
       return False


def runtest(test_name, test_func, successful, failed):
    if runtest1(test_func):
        successful.append(test_name)
    else:
        failed.append(test_name)

def runtests():
    printlevel = 1
    logging.logger.info(f"Executing {__name__} unit tests.")
    
    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    #
    tests = {"acquireMoleculeFiles": runAcquireMoleculeFilesTests, "acquireConfigurationFile": runAcquireConfigurationFileTests}

    #
    # Run
    #
    successful = []
    failed = []
    for test_name, test_func in tests.items():
        runtest(test_name, test_func, successful, failed)

    # 
    # Report
    #
    if not failed:        
        logging.logger.info(f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            logging.logger.error(f'Failed {test}.')

if __name__ == "__main__":
    runtests()