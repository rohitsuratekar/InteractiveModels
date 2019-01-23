from matplotlib import pylab as plt
from scipy.integrate import odeint

from analysis.helper import *


def test_para(row_no: int):
    enz = extract_enzyme_set("para.txt", row_no)
    ss_time_array = np.linspace(0, 2000, 5000)
    initial_condition = get_random_concentration(200)

    ss = odeint(get_equations(S_OPEN_2), initial_condition, ss_time_array,
                args=(enz, None))[-1]

    plc_base = enz[E_PLC].v
    for e in enz:
        if e != E_SOURCE:
            enz[e].k = enz[e].k / sum(ss)
            enz[e].v = enz[e].v / plc_base
        else:
            enz[e].k = enz[e].k / plc_base

    testing_time = np.linspace(0, 100, 5000)
    init = get_random_concentration(1)
    ss_new = odeint(get_equations(S_OPEN_2), init, testing_time,
                    args=(enz, None))

    for l in LIPID_NAMES:
        plt.plot(testing_time, ss_new[:, LIPID_NAMES.index(l)], label=l,
                 color=COLORS_PRIMARY[LIPID_NAMES.index(l)], linewidth=2)

    plt.legend(loc=0)
    plt.show()


def test_light(row_no: int):
    enz = extract_enzyme_set("para.txt", row_no)
    ss_time_array = np.linspace(0, 2000, 5000)
    initial_condition = get_random_concentration(200)

    ss = odeint(get_equations(S_OPEN_2), initial_condition, ss_time_array,
                args=(enz, None))[-1]

    plc_base = enz[E_PLC].v
    for e in enz:
        if e != E_SOURCE:
            enz[e].k = enz[e].k / sum(ss)
            enz[e].v = enz[e].v / plc_base
        else:
            enz[e].k = enz[e].k / plc_base

    testing_time = np.linspace(0, 100, 5000)
    init = get_random_concentration(1)
    ss_new = odeint(get_equations(S_OPEN_2), init, testing_time,
                    args=(enz, None))

    for l in LIPID_NAMES:
        plt.plot(testing_time, ss_new[:, LIPID_NAMES.index(l)], label=l,
                 color=COLORS_PRIMARY[LIPID_NAMES.index(l)], linewidth=2)

    plt.legend(loc=0)
    plt.show()
