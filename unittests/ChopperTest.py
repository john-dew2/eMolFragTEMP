from pathlib import Path

from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.chopper import Chopper

usr_dir = Path.cwd()
dblbond = usr_dir.joinpath("eMolFragTEMP/unittests/data/1098067-c2-c2-double-bond")
dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")
pdb = dataPath.joinpath("pbd")
mol = dataPath.joinpath("mol")

def runChop():
  pass

def runChopTests():
  pass
  #mol = utilities.getRDKitMolecule(mol2.joinpath("DB01060.mol2"))
  #linker, bricks = Chopper.chop(mol)
  #pass

def runProcessFragments():
  pass

def runProcessFragmentsTests():
  pass

def runChopAll():
  pass

def runChopAllTests():
  pass

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
    tests = {"chop": runChopTests, "chopall": runChopAllTests, "processFragments": runProcessFragmentsTests
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