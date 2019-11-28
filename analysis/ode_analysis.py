import numpy as np
from scipy.integrate import odeint

from analysis.helper import get_equations
from constants.namespace import *


class ODEAnalysis:
    def __init__(self, system: str, total_lipid: float, enzymes: dict,
                 feed_para: dict = None):
        self.system = system
        self.start_time = 0
        self.enz = enzymes
        self.end_time = 1000
        self.buffer_time = 2
        self.feed_para = feed_para
        self.recovery_time_array = None
        self.buffer_array = None
        self._time_array = None
        self._concentration_array = None
        self._steady_state = None
        self.stimulated_concentrations = None
        self.recovery_concentration = None
        self.total_lipid = total_lipid
        all_ratios = np.random.uniform(0, 1, 8)
        self.initial_condition = []
        for i in range(8):
            self.initial_condition.append(
                all_ratios[i] * self.total_lipid / sum(all_ratios))

    @property
    def time_array(self):
        return self._time_array

    @time_array.setter
    def time_array(self, array):
        if self.time_array is None:
            self._time_array = array
        else:
            self._time_array = np.concatenate([self._time_array,
                                               self._time_array[-1] + array])

    @property
    def concentration_array(self):
        return self._concentration_array

    @concentration_array.setter
    def concentration_array(self, array):
        if self.concentration_array is None:
            self._concentration_array = array
        else:
            self._concentration_array = np.concatenate((
                self._concentration_array, array), axis=0)

    @property
    def steady_state(self):
        return self._steady_state

    @steady_state.setter
    def steady_state(self, value):
        self._steady_state = value

    def normalize_enzymes(self, end_time: float = 2000):
        it = np.linspace(0, end_time, 5000)
        ss = odeint(get_equations(self.system), self.initial_condition, it,
                    args=(self.enz, None))[-1]
        plc_base = self.enz[E_PLC].v
        for e in self.enz:
            if e != E_SOURCE:
                self.enz[e].k = self.enz[e].k / sum(ss)
                self.enz[e].v = self.enz[e].v / plc_base
            else:
                self.enz[e].k = self.enz[e].k / plc_base

    def attain_steady_state(self):
        it = np.linspace(0, 2000, 5000)
        self.steady_state = \
            odeint(get_equations(self.system), self.initial_condition,
                   it, args=(self.enz, None))[-1]

    def initialize(self):
        it = np.linspace(0, 2000, 5000)
        self.attain_steady_state()
        if self.feed_para is not None:
            multi_factor = 1
            for key in self.feed_para:
                h = self.feed_para[key][F_HILL_COEFFICIENT]
                a = self.feed_para[key][F_MULTIPLICATION_FACTOR]
                c = self.feed_para[key][F_CARRYING_CAPACITY]
                feedback_substrate = self.steady_state[I_PIP2]
                t = self.feed_para[key][F_TYPE_OF_FEEDBACK]

                reg = 1 + pow((feedback_substrate / c), h)
                fed = 1 + a * pow((feedback_substrate / c), h)
                if t == FEEDBACK_POSITIVE:
                    multi_factor = reg / fed  # This is reversed
                elif t == FEEDBACK_NEGATIVE:
                    multi_factor = fed / reg  # This is reversed

                self.enz[key].v *= multi_factor

                # Adjust all steady states with new Vmax
                self.steady_state = \
                    odeint(get_equations(self.system), self.initial_condition,
                           it, args=(self.enz, self.feed_para))[-1]

        self.buffer_array = np.linspace(0, self.buffer_time,
                                        self.buffer_time * 10)

        self.time_array = self.buffer_array - self.buffer_array[-1]
        self.concentration_array = odeint(get_equations(self.system),
                                          self.steady_state,
                                          self.buffer_array,
                                          args=(self.enz, self.feed_para))

    def stimulate(self, depletion_percentage):
        sim_ss = [x for x in self.steady_state]
        amount = sim_ss[I_PIP2] * (100 - depletion_percentage) / 100
        sim_ss[I_DAG] = sim_ss[I_DAG] + sim_ss[I_PIP2] - amount
        sim_ss[I_PIP2] = amount
        self.stimulated_concentrations = sim_ss

    def scale_parameters_by(self, value: float):
        for e in self.enz:
            self.enz[e].v = value * self.enz[e].v

    def stimulate_pi4k(self, depletion_percentage):
        """
        This is hypothetical scenario. Use only for testing purpose
        """
        sim_ss = [x for x in self.stimulated_concentrations]
        amount = sim_ss[I_PI4P] * (100 - depletion_percentage) / 100
        sim_ss[I_DAG] = sim_ss[I_DAG] + sim_ss[I_PI4P] - amount
        sim_ss[I_PI4P] = amount
        self.stimulated_concentrations = sim_ss

    def recover(self, time: int):
        self.recovery_time_array = np.linspace(0, time, time * 100)
        self.time_array = self.recovery_time_array
        self.recovery_concentration = odeint(get_equations(self.system),
                                             self.stimulated_concentrations,
                                             self.recovery_time_array,
                                             args=(self.enz, self.feed_para))
        self.concentration_array = self.recovery_concentration
