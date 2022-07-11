
import sys
from Molecule import Molecule
from Options import Options
from pathlib import Path
from rdkit import Chem
import AcquireFiles
import AcquireMolecules

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    args = sys.argv
    initializer = Options()
    initializer.parseCommandLine(args)
    
    #Input System
    files = AcquireFiles.acquireFiles(initializer)
    dataset = AcquireMolecules.acquireMolecules(files)
    
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()