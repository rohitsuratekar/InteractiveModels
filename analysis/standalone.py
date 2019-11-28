import matplotlib.pylab as plt

from analysis.helper import *
from analysis.ode_analysis import ODEAnalysis
from SecretColors.palette import Palette

LIGHT_STIMULATION = 85
LIPIDS_TO_ANALYZE = [I_PIP2]

FEED_PARA = {
    # E_PIP5K: {
    #     F_TYPE_OF_FEEDBACK: FEEDBACK_NEGATIVE,
    #     F_HILL_COEFFICIENT: 1,
    #     F_MULTIPLICATION_FACTOR: 60,
    #     F_CARRYING_CAPACITY: 0.3,
    #     F_FEED_SUBSTRATE_INDEX: I_PIP2
    # }
    E_PI4K: {
        F_TYPE_OF_FEEDBACK: FEEDBACK_POSITIVE,
        F_HILL_COEFFICIENT: 1,
        F_MULTIPLICATION_FACTOR: 30,
        F_CARRYING_CAPACITY: 0.01,
        F_FEED_SUBSTRATE_INDEX: I_PIP2
    }
}


def provide_condition(base, scale):
    base.buffer_time = 0.5
    base.scale_parameters_by(scale)
    base.normalize_enzymes()
    base.initialize()
    base.stimulate(LIGHT_STIMULATION)
    # base.stimulate_pi4k(LIGHT_STIMULATION)  # Hypothetical scenario
    base.recover(5)


def plot_standalone():
    enz = extract_enzyme_set("para.txt")
    enz2 = extract_enzyme_set("para.txt", 0)  # Required because ODEAnalysis
    # enz2.get(E_PIP5K).v *= 0.1
    # enz2.get(E_PI4K).v *= 0.1
    # changes enzyme values
    base = ODEAnalysis(S_OPEN_2, 1, enz)
    feed = ODEAnalysis(S_OPEN_2, 1, enz2, FEED_PARA)

    provide_condition(base, 5)
    provide_condition(feed, 1)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p = Palette()

    for lipid in LIPIDS_TO_ANALYZE:
        scaled_array = base.concentration_array[:, lipid]
        scaled_array = scaled_array - min(scaled_array)
        scaled_array = [x * 100 / max(scaled_array) for x in scaled_array]
        ax.plot(base.time_array, scaled_array,
                label=LIPID_NAMES[lipid], color=COLORS_PRIMARY[lipid],
                )

        # scaled_feed = feed.concentration_array[:, lipid]
        # scaled_feed = [x * 100 / max(scaled_feed) for x in scaled_feed]
        # ax.plot(feed.time_array, scaled_feed,
        #         label=LIPID_NAMES[lipid] + "(with Feedback)",
        #         color=COLORS_PRIMARY[lipid])

    ax.set_ylabel("Scaled Concentration")
    ax.set_xlabel("Time (arbitrary)")
    ax.legend(loc=0)
    plt.axvline(0, linestyle="--", color="k")
    plt.axhline(0, linestyle="--", color="k")
    plt.show()


def plot_test():
    enz = extract_enzyme_set("para.txt", 0)
    enz2 = extract_enzyme_set("para.txt", 0)  # Required because ODEAnalysis
    # changes enzyme values
    base = ODEAnalysis(S_OPEN_2, 1, enz)
    # feed = ODEAnalysis(S_OPEN_2, 1, enz2, FEED_PARA)

    provide_condition(base, 6)
    # provide_condition(feed)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for lipid in LIPIDS_TO_ANALYZE:
        scaled_array = base.concentration_array[:, lipid]
        scaled_array = [x * 100 / max(scaled_array) for x in scaled_array]
        ax.plot(base.time_array, scaled_array,
                label=LIPID_NAMES[lipid], color=COLORS_PRIMARY[lipid],
                )
        # ax.plot(feed.time_array, feed.concentration_array[:, lipid],
        #         label=LIPID_NAMES[lipid] + "(with Feedback)",
        #         color=COLORS_PRIMARY[lipid])

    ax.set_ylabel("Scaled Concentration")
    ax.set_xlabel("Time (arbitrary)")
    ax.legend(loc=0)
    plt.axvline(0, linestyle="--", color="k")
    plt.show()


def run():
    plot_standalone()
