from pathlib import Path

from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.chopper import Connectivity, Deconstructor
from eMolFragTEMP.src.input import AcquireMolecules, AcquireFiles, Options, Configuration

usr_dir = Path.cwd()
dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")
pdb = dataPath.joinpath("pbd")
mol = dataPath.joinpath("mol")

def getListofFiles(arguments):
    arg = utilities.createParser(arguments)
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arg)
    return AcquireFiles.acquireMoleculeFiles(initializer)

def runCompute(usr_file, expec):
  mol = AcquireMolecules.acquireMolecules(getListofFiles(f"-i {usr_file} -o output/".split(" ")))
  mol = mol[0].rdkitObject
  print(mol)
  snips = Deconstructor.deconstruct(mol)
  
  snips = snips[2]
  #print(snips)
  Connectivity.compute(mol, snips)
  assert Connectivity.compute(mol, snips) == expec

def runComputeTests():
  runCompute(mol2, 0)

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
    
def runtests(printlevel):

    utilities.emit(printlevel, f'Executing {__name__} unit tests.')

    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    #
    tests = {"compute": runComputeTests
            }

    #
    # Run
    #
    successful = []
    failed = []
    for (test_name, test_func) in tests.items():
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
    runtests(0)