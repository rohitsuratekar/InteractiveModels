import json

from constants.namespace import *
from models.biology import Enzyme
from models.systems import open2


def extract_enzyme_set(filename: str, row: int = 0) -> dict:
    """
    Only extract first enzyme value.
    You can provide option option of providing row number for line
    :param filename:
    :param row: Row number of parameter value
    :return: Dictionary of Enzymes
    """
    count = 0
    with open(filename) as f:
        for line in f:
            if count == row:
                raw_enz = json.loads(line.split(":", 1)[1])["Enzymes"]
                enz = {}
                for key in raw_enz:
                    enz[key] = Enzyme.make(key, raw_enz[key])
                return enz
            else:
                count += 1


def get_equations(system: str):
    """
    Returns equation of specific system
    :param system: topology or known model
    :return: set of equation function
    """
    if system == S_OPEN_2:
        return open2
    else:
        raise Exception("No such system found :%s" % system)


COLORS_PRIMARY = ["#00a78f", "#db7c00", "#12a3b4", "#ff509e", "#5392ff",
                  "#00aa5e", "#b3901f", "#b07ce8"]

LIPID_NAMES = [L_PMPI, L_PI4P, L_PIP2, L_DAG, L_PMPA, L_ERPA, L_CDPDAG, L_ERPI]
