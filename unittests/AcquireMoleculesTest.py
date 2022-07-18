import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.input import AcquireMolecules, AcquireFiles, Configuration, Options

usr_dir = Path.cwd()
dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")
pdb = dataPath.joinpath("pbd")
mol = dataPath.joinpath("mol")

def runAcquireMoleculesTests():

    testPaths = [mol2, smi, pdb, mol]
    #Tests if mol2 and smi are taken
    for filePath in testPaths:
        runAcquireMolecules(getListofFiles(f"-i {filePath} -o output/".split(" ")), 5)
        
    #Will recognize the files as bad
    runAcquireMolecules(getListofFiles(f"-i {sdf} -o output/".split(" ")), 0)

def getListofFiles(arguments):
    print(arguments)
    arg = utilities.createParser(arguments)
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arg)
    return AcquireFiles.acquireMoleculeFiles(initializer)

def runAcquireMolecules(files, expec):
    retFiles = AcquireMolecules.acquireMolecules(files)
    assert len(retFiles) == expec

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
    utilities.emit(printlevel, f"Executing {__name__} unit tests.")
    
    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    #
    tests = {"acquireMolecules": runAcquireMoleculesTests}

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
        utilities.emit(printlevel, f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            utilities.emit(printlevel+1, f'Failed {test}.')

if __name__ == "__main__":
    runtests()