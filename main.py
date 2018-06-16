from analysis.interactive import interact
from analysis.interactive_mutant import interact as mt_interact
from analysis.standalone import plot_standalone


def stand():
    plot_standalone()


def inte():
    interact()


def mut():
    mt_interact()


if __name__ == "__main__":
    interact()
