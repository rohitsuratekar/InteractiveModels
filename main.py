from analysis.interactive_light import interact
from analysis.interactive_mutant import interact as mt_interact
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


if __name__ == "__main__":
    inter()
