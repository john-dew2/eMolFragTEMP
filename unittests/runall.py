#
# A 'main' function for all unit tests of emolfrag
#
# Tests are run in increasing levels of project dependence
#     i.e., files with little intra-project dependent are run first
#

import MoleculeDatabaseTest
import ChopperTest

def runtests():

    emit(0, f'Executing eMolFrag v2.0 unit tests.')
    
    MoleculeDatabaseTest.runtests(1)

    ChopperTest.runtests()

        
if __name__ == "__main__":
    runtests()