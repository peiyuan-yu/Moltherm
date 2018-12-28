from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor

from moltherm.compute.processing import MolThermDataProcessor
from atomate.qchem.database import QChemCalcDb


__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "December 2018"


class ReactionRanker:
    """
    This class includes several ranking algorithms which can distinguish between
    different reactions. For each algorithm, flexibility in parameters is
    allowed, and in general, the same set of parameters can be used for each
    algorithm.

    At present, the allowed set of parameters includes the following:
    bp - normal boiling point
    mp - normal melting point
    vp - vapor pressure
    solubility - aqueous solubility
    low_kow - octanol-water partition coefficient
    """

    def __init__(self):
        pass

    def pareto_ranking(self):
        pass

    def heuristic_ranking(self):
        pass