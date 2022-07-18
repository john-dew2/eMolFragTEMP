import sys
from pathlib import Path
from eMolFragTEMP.unittests import utilities

#from ..src.input import Configuration, Options
from eMolFragTEMP.src.input import Configuration, Options

usr_dir = Path.cwd()
config_files = usr_dir.joinpath("eMolFrag2/unittests/data/config-files")

dataPath = usr_dir.joinpath("eMolFrag2/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")

def createParser(arguments):
  if (arguments == None):
    return None

  from argparse import ArgumentParser
  parser = ArgumentParser(description='eMolFrag2')
  parser.add_argument("-i",
  type=str,
  help='Set the input path')
  
  parser.add_argument("-o",
  type=str,
  help='Set the output path')
  
  #parser.add_argument("-m",
  #type=int, choices=range(0,3),
  #help='Set the execution type')
  
  parser.add_argument("-c",
  type=int, choices=range(0,3),
  help='Set the output type')
  
  
  #print("pass")
  args = parser.parse_args(arguments)
  print(args)
  return parser.parse_args(arguments)

def runReadCommandLine(arguments, expec_result):
    arg = createParser(arguments)
    initializer = Options.Options()
    initializer = Configuration.readCommandLine(initializer, arg)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadCommandLineTests():
    print("yuh")
    #Normal
    #runReadCommandLine(f"-i {mol2} -o output/ -c 1".split(" "), True)
    
    #Empty
    #runReadCommandLine(None, False)
    
    
def runReadConfigurationFile(config_file, expec_result):
  assert Configuration.readConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    print("yuh")
    #Files with a comment
    runReadConfigurationFile(config_files.joinpath("comment1.txt"), ['\n!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(config_files.joinpath("comment2.txt"), ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1', ''])
    runReadConfigurationFile(config_files.joinpath("comment3.txt"), ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1\n\n'])
    runReadConfigurationFile(config_files.joinpath("comment4.txt"), [''])

    #Normal Files
    runReadConfigurationFile(config_files.joinpath("normal1.txt"),  ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(config_files.joinpath("normal2.txt"),  ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '0', '-c', '0'])

    #Empty File
    runReadConfigurationFile(config_files.joinpath("empty1.txt"), [])

def runReadConfigurationInput(arguments, expec_result):
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arguments)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadConfigurationInputTests():
    print("yuh")
    #Normal processing
    #runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {mol2} -o output/ -p 16 -c 1".split(" "), True)
    #runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {smi} -o output/".split(" "), True)
    
    #No input or output
    #runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -o output/ -p 16 -c 1".split(" "), False)
    #runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {smi} -p 16 -c 1".split(" "), False)
    #runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -p 16 -c 1".split(" "), False)
    
    #Not enough arguments 
    #runReadConfigurationInput(f"!python".split(" "), False)
    
    #Files only
        #Files only
    #runReadConfigurationInput(f"-m {config_files.joinpath('comment1.txt')}".split(" "), True)
    #runReadConfigurationInput(f"-m {config_files.joinpath('comment2.txt')}".split(" "), True)
    #runReadConfigurationInput(f"-m {config_files.joinpath('comment3.txt')}".split(" "), True)
    #runReadConfigurationInput(f"-m {config_files.joinpath('comment4.txt')}".split(" "), False)
    #runReadConfigurationInput(f"-m {config_files.joinpath('empty1.txt')}".split(" "), False)
    #runReadConfigurationInput(f"-m {config_files.joinpath('normal1.txt')}".split(" "), True)
    #runReadConfigurationInput(f"-m {config_files.joinpath('normal2.txt')}".split(" "), True)
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