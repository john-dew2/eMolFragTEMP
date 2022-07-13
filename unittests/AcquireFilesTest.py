import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.input import AcquireFiles, Configuration, Options

usr_dir = Path.cwd()
dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")

def runAcquireMoleculeFilesTests():
    
    testPaths = [mol2, smi]
    for filePath in testPaths:
        runAcquireMoleculeFiles(f"-m eMolFragTEMP.src.eMolFrag -i {filePath} -o output/".split(" "), 5)
    runAcquireMoleculeFiles(f"-m eMolFragTEMP.src.eMolFrag -i {sdf} -o output/".split(" "), 0)
    runAcquireMoleculeFiles(f"-m eMolFragTEMP.src.eMolFrag -i directory/doesnt/exist -o output/".split(" "), 0)
        
def runAcquireMoleculeFiles(arguments, expec):
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arguments)
    files = AcquireFiles.acquireMoleculeFiles(initializer)
    if (files == None): files = []
    assert len(files) == expec
    
def runAcquireConfigurationFileTests():
    runAcquireConfigurationFile(mol2.joinpath("DB00607.mol2"), mol2.joinpath("DB00607.mol2"))
    runAcquireConfigurationFile(smi.joinpath("DB00607.smi"), smi.joinpath("DB00607.smi"))
    runAcquireConfigurationFile(sdf.joinpath("DB00607.sdf"), sdf.joinpath("DB00607.sdf"))
    runAcquireConfigurationFile(mol2.joinpath("DB12345.mol2"), None)

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
    utilities.emit(printlevel, f"Executing {__name__} unit tests.")
    
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
        utilities.emit(printlevel, f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            utilities.emit(printlevel+1, f'Failed {test}.')

if __name__ == "__main__":
    runtests()