from analysis.interactive_light import interact
from analysis.interactive_mutant import interact as mt_interact
from analysis.pi_analysis import plot_pi
from analysis.standalone import plot_standalone
from analysis.test import test_para


def stand():
    plot_standalone()


def inter():
    interact()


def mut():
    mt_interact()


def test(row):
    test_para(row)


def pi():
    plot_pi()


if __name__ == "__main__":
    mut()
