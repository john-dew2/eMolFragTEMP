import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.utilities import logging
#from ..src.input import Configuration, Options
from eMolFragTEMP.src.input import Configuration, Options

usr_dir = Path.cwd()
config_files = usr_dir.joinpath("eMolFragTEMP/unittests/data/config-files")

dataPath = usr_dir.joinpath("eMolFragTEMP/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")

def runReadCommandLine(arguments, expec_result):
    arg = utilities.createParser(arguments)
    initializer = Options.Options()
    initializer = Configuration.readCommandLine(initializer, arg)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadCommandLineTests():
    #Normal
    runReadCommandLine(f"-i {mol2} -o output/".split(" "), True)
    
    #Empty
    runReadCommandLine(None, False)

     
def runReadConfigurationFile(config_file, expec_result):
  assert Configuration.readConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    #Files with a comment
    runReadConfigurationFile(config_files.joinpath("comment1.txt"), ['-i', '/content/eMolFrag2/unittests/data/db-files', '-o', 'output/', '-u', '-indiv'])
    runReadConfigurationFile(config_files.joinpath("comment2.txt"), ['-i', '/content/eMolFrag2/unittests/data/db-files', '-o', 'output/', '-u', '-indiv'])
    #runReadConfigurationFile(config_files.joinpath("comment3.txt"), ['-i', '/content/eMolFrag2/unittests/data/db-files', '-o', 'output/', '-u', '-indiv\n'])
    runReadConfigurationFile(config_files.joinpath("comment4.txt"), [''])
    #Normal Files
    runReadConfigurationFile(config_files.joinpath("normal1.txt"),  ['-i', '/content/eMolFrag2/unittests/data/db-files', '-o', 'output/', '-u', '-indiv'])
    runReadConfigurationFile(config_files.joinpath("normal2.txt"),  ['-i', '/content/eMolFrag2/unittests/data/db-files', '-o', 'output/'])

    #Empty File
    runReadConfigurationFile(config_files.joinpath("empty1.txt"), [])

def runReadConfigurationInput(arguments, expec_result):
    arg = utilities.createParser(arguments)
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arg)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadConfigurationInputTests():
    #Normal processing
    runReadConfigurationInput(f"-i {mol2} -o output/".split(" "), True)
    runReadConfigurationInput(f"-i {smi} -o output/".split(" "), True)
    #No input or output
    runReadConfigurationInput(f"-o output/".split(" "), False)
    runReadConfigurationInput(f"-i {smi}".split(" "), False)
    runReadConfigurationInput(f"-indiv".split(" "), False)
    #Not enough arguments 
    #runReadConfigurationInput(f"".split(" "), False)
    
    #Files only
    runReadConfigurationInput(f"-i {config_files.joinpath('comment1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-i {config_files.joinpath('comment2.txt')}".split(" "), True)
    runReadConfigurationInput(f"-i {config_files.joinpath('comment3.txt')}".split(" "), True)
    runReadConfigurationInput(f"-i {config_files.joinpath('comment4.txt')}".split(" "), False)
    runReadConfigurationInput(f"-i {config_files.joinpath('empty1.txt')}".split(" "), False)
    runReadConfigurationInput(f"-i {config_files.joinpath('normal1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-i {config_files.joinpath('normal2.txt')}".split(" "), True)
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
        logging.logger.info(f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            logging.logger.error(f'Failed {test}.')

if __name__ == "__main__":
    runtests()