import matplotlib.pylab as plt

from analysis.helper import *
from analysis.ode_analysis import ODEAnalysis

MUTANT = E_SINK
PARAMETER_SET = 2


def plot_pi():
    enz = extract_enzyme_set("para.txt", PARAMETER_SET)

    # changes enzyme values
    base = ODEAnalysis(S_OPEN_2, 1, enz)

    if PARAMETER_SET == 0:
        base.normalize_enzymes()

    base.attain_steady_state()
    s = [x for x in base.steady_state]
    s.append(
        base.steady_state[I_PMPA] + base.steady_state[I_ERPA])
    s.append(
        base.steady_state[I_PMPI] + base.steady_state[I_ERPI])

    base.enz[MUTANT].v *= 0.1
    if MUTANT == E_SOURCE:
        base.enz[MUTANT].k *= 0.1

    base.attain_steady_state()
    mt = [x for x in base.steady_state]
    mt.append(base.steady_state[I_PMPA] + base.steady_state[I_ERPA])
    mt.append(base.steady_state[I_PMPI] + base.steady_state[I_ERPI])

    figure = plt.figure()
    ax = figure.add_subplot(111)
    ind = np.arange(len(LIPID_NAMES) + 2)  # the x locations for the groups
    width = 0.35  # the width of the bars

    ax.bar(ind, np.asanyarray(mt) / np.asanyarray(s), width, color="#00b6cb")
    ax.set_xticks(ind)
    k = [round(x, 3) for x in np.asanyarray(mt) / np.asanyarray(s)]
    km = ""
    for s in k:
        km += str(s) + "\t"
    print(km)

    plt.axhline(1, linestyle="--", color="k")
    # plt.yscale("log")
    nm = [L_PMPI, L_PI4P, L_PIP2, L_DAG, L_PMPA, L_ERPA, L_CDPDAG,
          L_ERPI, "PA$_{total}$", "PI$_{total}$"]
    ax.set_xticklabels(nm)
    ax.set_ylim(0, 2)
    ax.set_ylabel("$\\frac{Mutant}{Wild Type}$ steady states")

    # plt.show()
