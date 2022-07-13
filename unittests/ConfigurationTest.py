import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities

#from ..src.input import Configuration, Options
from eMolFragTEMP.src.input import Configuration, Options

usr_dir = Path.cwd()

def runReadCommandLine(arguments, expec_result):
    initializer = Options.Options()
    initializer = Configuration.readCommandLine(initializer, arguments)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadCommandLineTests():
    runReadCommandLine("-m eMolFragTEMP.src.eMolFrag -i {usr_dir.joinpath('eMolFragTEMP/test/mol2-test')} -o output/ -p 16 -c 1".split(" "), True)
    runReadCommandLine(None, False)
    
    
def runReadConfigurationFile(config_file, expec_result):
    assert Configuration.readConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/comment1.txt"), ['\n!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/comment2.txt"), ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1', ''])
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/comment3.txt"), ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1\n\n'])
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/comment4.txt"), [''])
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/normal1.txt"),  ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files/normal2.txt"),  ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '0', '-c', '0'])

def runReadConfigurationInput(arguments, expec_result):
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arguments)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadConfigurationInputTests():
    runReadConfigurationInput(f"-m eMolFragTEMP.src.eMolFrag -i {usr_dir.joinpath('eMolFragTEMP/test/mol2-test')} -o output/ -p 16 -c 1".split(" "), True)
    runReadConfigurationInput(f"-m eMolFragTEMP.src.eMolFrag -o output/ -p 16 -c 1".split(" "), False)
    runReadConfigurationInput(f"-m eMolFragTEMP.src.eMolFrag -i {usr_dir.joinpath('eMolFragTEMP/test/mol2-test')} -p 16 -c 1".split(" "), False)
    runReadConfigurationInput(f"-m eMolFragTEMP.src.eMolFrag -i {usr_dir.joinpath('eMolFragTEMP/test/mol2-test')} -o output/".split(" "), True)
    runReadConfigurationInput(f"!python".split(" "), False)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/comment1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/comment2.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/comment3.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/comment4.txt')}".split(" "), False)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/normal1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {usr_dir.joinpath('eMolFragTEMP/unittests/data/config-files/normal2.txt')}".split(" "), True)

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