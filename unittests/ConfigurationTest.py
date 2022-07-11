import sys
from pathlib import Path
import utilities


def runReadCommandLine(arguments, expec_result):
    initializer = Options()
    assert Configuration.readCommandLine(initializer, arguments) == expec_result

def runReadCommandLineTests():
    #runReadCommandLine(argument, expec_result)
    pass
    
def runReadConfigurationFile(config_file, expec_result):
    assert Configuration.ReadConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    runReadConfigurationFile("./data/config-files/comment1.txt", ['eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("./data/config-files/comment2.txt", ['eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("./data/config-files/comment3.txt", ['eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("./data/config-files/comment4.txt", [])
    runReadConfigurationFile("./data/config-files/normal1.txt",  ['eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("./data/config-files/normal2.txt",  ['eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '0', '-c', '0'])

def runReadConfigurationInput(arguments, expec_result):
    initializer = Options()
    assert Configuration.ReadConfigurationInput(initializer, arguments) == expec_result

def runReadConfigurationInputTests():
    pass


#
#
# Unit test initiation functionality
#
#
def runtest(func):

    try:
       func()
       return True

    except:
       return False


def runtest(test_name, test_func, successful, failed):
    if runtest(test_func):
        successful.append(test_name)
    else:
        failed.append(test_name)

def runtests():
    printlevel = 1
    utilities.emit(f"Executing {__name__} unit tests.")
    
    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    #
    tests = {"readCommandLine": runReadCommandLineTests, "readConfigurationFile": runReadConfigurationFileTests, "readConfigurationInput": runReadConfigurationInputTests
            }

    #
    # Run
    #
    successful = []
    failed = []
    for (test_name, test_func) in tests:
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