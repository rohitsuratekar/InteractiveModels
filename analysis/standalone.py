import matplotlib.pylab as plt

from analysis.helper import *
from analysis.ode_analysis import ODEAnalysis

LIGHT_STIMULATION = 85
LIPIDS_TO_ANALYZE = [I_PIP2, I_PI4P]

FEED_PARA = {
    E_PIP5K: {
        F_TYPE_OF_FEEDBACK: FEEDBACK_NEGATIVE,
        F_HILL_COEFFICIENT: 2,
        F_MULTIPLICATION_FACTOR: 8,
        F_CARRYING_CAPACITY: 0.1,
        F_FEED_SUBSTRATE_INDEX: I_PIP2
    }
}


def provide_condition(base):
    base.buffer_time = 0.5
    base.normalize_enzymes()
    base.initialize()
    base.stimulate(LIGHT_STIMULATION)
    base.recover(2)


def plot_standalone():
    enz = extract_enzyme_set("para.txt")
    enz2 = extract_enzyme_set("para.txt")  # Required because ODEAnalysis
    # changes enzyme values
    base = ODEAnalysis(S_OPEN_2, 1, enz)
    feed = ODEAnalysis(S_OPEN_2, 1, enz2, FEED_PARA)

    provide_condition(base)
    provide_condition(feed)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for lipid in LIPIDS_TO_ANALYZE:
        ax.plot(base.time_array, base.concentration_array[:, lipid],
                label=LIPID_NAMES[lipid], color=COLORS_PRIMARY[lipid],
                linestyle="--")
        ax.plot(feed.time_array, feed.concentration_array[:, lipid],
                label=LIPID_NAMES[lipid] + "(with Feedback)",
                color=COLORS_PRIMARY[lipid])

    ax.set_ylabel("Scaled Concentration")
    ax.set_xlabel("Time (arbitrary)")
    ax.legend(loc=0)
    plt.axvline(0, linestyle="--", color="k")
    plt.show()
