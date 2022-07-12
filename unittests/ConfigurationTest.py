import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities

#from ..src.input import Configuration, Options
from eMolFragTEMP.src.input import Configuration, Options


def runReadCommandLine(arguments, expec_result):
    initializer = Options()
    assert Configuration.readCommandLine(initializer, arguments) == expec_result

def runReadCommandLineTests():
    #runReadCommandLine(argument, expec_result)
    pass
    
def runReadConfigurationFile(config_file, expec_result):
    assert Configuration.readConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/comment1.txt", ['\n!python', 'eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/comment2.txt", ['!python', 'eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1', ''])
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/comment3.txt", ['!python', 'eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1\n\n'])
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/comment4.txt", [''])
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/normal1.txt",  ['!python', 'eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile("/content/eMolFragTEMP/unittests/data/config-files/normal2.txt",  ['!python', 'eMolFrag-2/src/eMolFrag.py', '-i', '/content/input', '-o', 'output/', '-p', '0', '-c', '0'])

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
    tests = {"readCommandLine": runReadCommandLineTests, "readConfigurationFile": runReadConfigurationFileTests, "readConfigurationInput": runReadConfigurationInputTests}

    #
    # Run
    #
    successful = []
    failed = []
    for test_name, test_func in tests.items():
        print(f"Running {test_name}...")
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