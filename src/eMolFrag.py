
import sys
from eMolFragTEMP.src.representation import Molecule
from pathlib import Path
from rdkit import Chem
from eMolFragTEMP.src.input import AcquireFiles, AcquireMolecules, Configuration, Options
from eMolFragTEMP.unittests import utilities
from eMolFragTEMP.src.utilities import logging
from eMolFragTEMP.src.chopper import Chopper

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    #args = sys.argv
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, ARGS)

    # Get files
    files = AcquireFiles.acquireMoleculeFiles(initializer)    
    logging.logger.info(f'{len(files)} files to be processed.')
    
    # Get molecules
    molecules = AcquireMolecules.acquireMolecules(files)
    logging.logger.info(f'{len(molecules)} molecules to be chopped.')
 
    # CHOP
    brick_db, linker_db = Chopper.chopall(molecules)
    
    # Output fragments
    logging.logger.info(f'{brick_db.numUnique()} unique bricks among {brick_db.numAllMolecules()} bricks')
    logging.logger.info(f'{linker_db.numUnique()} unique linkers among {linker_db.numAllMolecules()} linkers')

    print(brick_db)    
    print(linker_db)
	
	





if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='eMolFrag2')
  parser.add_argument("-i",
  type=str,
  help='Set the input path')
  
  parser.add_argument("-o",
  type=str,
  help='Set the output path')
  
  parser.add_argument("-u",
  type=int, choices = True, False
  help='Set the execution type')
  
  parser.add_argument("-indiv",
  type=int, choices = True, False
  help='Set the output type')
  
  
  

  ARGS = parser.parse_args()
  
  if ((ARGS.i == None) or (ARGS.o == None)):
    utilities.emit(0, f"Every command must iclude '-i' and '-o'")
  else:
  #logging.basicConfig(level=args.loglevel)
    main()