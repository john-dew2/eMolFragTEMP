import sys
from pathlib import Path
from eMolFrag2.unittests import utilities
from eMolFrag2.src.input import Configuration, Options

def runAcquireMoleculeFilesTests():

def runAcquireMoleculeFiles():

def runAcquireConfigurationFileTests():

def runAcquireConfigurationFile():






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