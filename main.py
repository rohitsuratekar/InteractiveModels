from analysis.interactive_light import interact
from analysis.interactive_mutant import interact as mt_interact
from analysis.pi_analysis import plot_pi
from analysis.standalone import run
from analysis.test import test_light


def stand():
    run()


def inter():
    interact()


def mut():
    mt_interact()


def test(row):
    test_light(row)


def pi():
    plot_pi()


if __name__ == "__main__":
    stand()
