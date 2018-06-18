import sys
from functools import partial

import matplotlib.pylab as plt
from PyQt5.QtWidgets import QApplication, QDialog, QGridLayout, QComboBox, \
    QLabel, QHBoxLayout, QVBoxLayout, QGroupBox, QDoubleSpinBox, QPushButton
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar

from analysis.helper import *
from analysis.ode_analysis import ODEAnalysis

E_NONE = "None"
MUTANT_DEPLETION = 0.1
PARAMETER_NO = 1


class MySpinBox(QDoubleSpinBox):
    def __init__(self, *args):
        QDoubleSpinBox.__init__(self, *args)
        self.setDecimals(4)
        self.setRange(0.0001, 10)
        self.setSingleStep(0.05)
        self.setRange(0.0001, 1000)

    def textFromValue(self, value):
        return "%.4f" % value


def get_ratios(wt, mt):
    wt_pi = wt[I_PMPI] + wt[I_ERPI]

    mt_pi = mt[I_PMPI] + mt[I_ERPI]

    mt_pa = mt[I_PMPA] + mt[I_ERPA]
    wt_pa = wt[I_PMPA] + wt[I_ERPA]

    return (mt[I_DAG] / mt_pi) / (wt[I_DAG] / wt_pi), (mt_pa / mt_pi) / (
            wt_pa / wt_pi)


class MutantWindow(QDialog):
    def __init__(self, parent=None):
        super(MutantWindow, self).__init__(parent)
        self.setWindowTitle('Mutant Interactive Analysis')
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.show_standard = False
        self.para_dict = {}
        self.mutant_enzyme = E_NONE
        self.enz_box = None

        enz = extract_enzyme_set("para.txt", PARAMETER_NO)
        self.base = ODEAnalysis(S_OPEN_2, 1, enz)
        if PARAMETER_NO == 0:
            self.base.normalize_enzymes()
        self.base.attain_steady_state()

        self.parameter_box = QGridLayout()
        self.ss_label = QLabel("Steady State Ratios: PIP2: 0.05, PI4P:0.05, "
                               "DAG: 0.008, PA: 0.1677")

        self.grid = QGridLayout()
        self.plot()
        self.set_layout()

    @property
    def enzyme_parameters(self):
        p_box = QGridLayout()
        col = 0
        row = 1
        for i in ENZYME_NAMES:
            main_box_layout = QGroupBox(i.upper())
            main_box = QVBoxLayout()
            main_box_layout.setLayout(main_box)
            main_box_layout.setContentsMargins(1, 10, 10, 1)

            box = QHBoxLayout()
            box.addWidget(QLabel("Vmax"))
            v = MySpinBox()
            v.setValue(self.base.enz[i].v)
            if i == E_PLC or i == E_SOURCE:
                v.setEnabled(False)
            box.addWidget(v)
            box.addWidget(QLabel("Km"))
            k = MySpinBox()
            k.setValue(self.base.enz[i].k)
            v.valueChanged.connect(partial(self.on_parameter_change, i, "v",
                                           v))
            k.valueChanged.connect(partial(self.on_parameter_change, i, "k",
                                           k))
            self.para_dict[i] = {
                "v": v, "k": k
            }
            box.addWidget(k)
            box.addStretch()
            if col == 3:
                row += 1
                col = 1
            else:
                col += 1
            main_box.addLayout(box)

            p_box.addWidget(main_box_layout, row, col)
        self.parameter_box = p_box
        return p_box

    @property
    def mutant_list(self):
        e_box = QHBoxLayout()
        e_box.setContentsMargins(30, 1, 1, 1)
        e_box.addStretch()
        e_box.addWidget(QLabel("Mutant Enzyme : "))
        self.enz_box = QComboBox(self)
        self.enz_box.addItem(E_NONE)
        for e in ENZYME_NAMES:
            self.enz_box.addItem(e)
            self.enz_box.activated[str].connect(self.on_mutant_selection)
        e_box.addWidget(self.enz_box)
        e_box.addStretch()
        return e_box

    @property
    def reset_btn(self):
        btn = QPushButton("Reset")
        btn.clicked[bool].connect(self.on_reset_click)
        return btn

    @property
    def print_btn(self):
        btn = QPushButton("Print Parameters")
        btn.clicked[bool].connect(self.print)
        return btn

    @property
    def toggle_plot(self):
        btn = QPushButton("Toggle Plot")
        btn.clicked[bool].connect(self.on_toggle_plot)
        return btn

    def on_toggle_plot(self):
        if self.show_standard:
            self.show_standard = False
        else:
            self.show_standard = True
        self.plot()

    def print(self):
        enz = {}
        for k in self.base.enz:
            enz[k] = self.base.enz[k].properties

        data = {"Enzymes": enz}
        print(json.dumps(data))

    def reset(self):
        enz = extract_enzyme_set("para.txt", PARAMETER_NO)
        self.base = ODEAnalysis(S_OPEN_2, 1, enz)
        if PARAMETER_NO == 0:
            self.base.normalize_enzymes()
        self.base.attain_steady_state()
        for key in self.para_dict:
            self.para_dict[key]["v"].setValue(self.base.enz[key].v)
            self.para_dict[key]["k"].setValue(self.base.enz[key].k)

    def on_reset_click(self):
        self.mutant_enzyme = E_NONE
        self.enz_box.setCurrentIndex(0)
        self.reset()
        self.plot()

    def on_parameter_change(self, enzyme, pro, new_val):
        if pro == "v":
            self.base.enz[enzyme].v = new_val.value()
        elif pro == "k":
            self.base.enz[enzyme].k = new_val.value()
        self.plot()

    def on_mutant_selection(self, value):
        self.mutant_enzyme = value
        self.reset()
        self.plot()

    def plot(self):
        self.figure.clear()
        self.base.attain_steady_state()
        s = self.base.steady_state
        p = s[I_PMPI] + s[I_ERPI]
        e_pip2 = abs((s[I_PIP2] / p) / 0.05 - 1)
        e_pi4p = abs((s[I_PI4P] / p) / 0.05 - 1)
        e_dag = abs((s[I_DAG] / p) / 0.008 - 1)
        e_pa = abs(((s[I_PMPA] + s[I_ERPA]) / p) / 0.1677 - 1)
        e_cdpdag = abs((s[I_CDPDAG] / p) / 0.001 - 1)
        error = e_pip2 + e_pi4p + e_dag + e_pa + e_cdpdag
        self.ss_label.setText("Steady State Ratios \n(w.r.t. Total PI): "
                              "\n\nPIP2: %.4f, \nPI4P: %.4f, "
                              "\nDAG: %.4f, \nPA: %.4f, \nCDPDAG: "
                              "%.4f\n\n\nTotal WT Error: "
                              "%.3f %% \n\nPMPI/ERPI: "
                              "%.4f\nPMPA/ERPA: %.4f" %
                              (s[I_PIP2] / p, s[I_PI4P] / p, s[I_DAG] / p,
                               (s[I_PMPA] + s[I_ERPA]) / p, s[I_CDPDAG] / p,
                               error * 100,
                               s[I_PMPI] / s[I_ERPI], s[I_PMPA] / s[I_ERPA]))

        if self.show_standard:
            self.standard_plot()
        else:
            self.individual_plot()
        self.canvas.draw()

    def individual_plot(self):
        s = [x for x in self.base.steady_state]
        s.append(
            self.base.steady_state[I_PMPA] + self.base.steady_state[I_ERPA])
        s.append(
            self.base.steady_state[I_PMPI] + self.base.steady_state[I_ERPI])

        if self.mutant_enzyme != E_NONE:
            self.base.enz[self.mutant_enzyme].v *= MUTANT_DEPLETION

        self.base.attain_steady_state()
        mt = [x for x in self.base.steady_state]
        mt.append(
            self.base.steady_state[I_PMPA] + self.base.steady_state[I_ERPA])
        mt.append(
            self.base.steady_state[I_PMPI] + self.base.steady_state[I_ERPI])

        ax = self.figure.add_subplot(111)
        ind = np.arange(len(LIPID_NAMES) + 2)  # the x locations for the groups
        width = 0.35  # the width of the bars

        ax.bar(ind, np.asanyarray(mt) / np.asanyarray(s), width)
        ax.set_xticks(ind)

        plt.axhline(1, linestyle="--")
        # plt.yscale("log")
        nm = [L_PMPI, L_PI4P, L_PIP2, L_DAG, L_PMPA, L_ERPA, L_CDPDAG,
              L_ERPI, "PA", "PI"]
        ax.set_xticklabels(nm)

        if self.mutant_enzyme != E_NONE:
            self.base.enz[self.mutant_enzyme].v /= MUTANT_DEPLETION

    def standard_plot(self):
        s = [x for x in self.base.steady_state]
        self.base.enz[E_DAGK].v *= MUTANT_DEPLETION

        self.base.attain_steady_state()
        mt_rdga = [x for x in self.base.steady_state]

        self.base.enz[E_DAGK].v /= MUTANT_DEPLETION

        self.base.enz[E_LAZA].v *= MUTANT_DEPLETION

        self.base.attain_steady_state()
        mt_laza = [x for x in self.base.steady_state]

        self.base.enz[E_LAZA].v /= MUTANT_DEPLETION

        ax = self.figure.add_subplot(111)
        ind = np.arange(4)  # the x locations for the groups
        width = 0.35  # the width of the bars

        exp_observations = [1, 1, 1, 2.5]

        vals1 = get_ratios(s, mt_rdga)
        vals2 = get_ratios(s, mt_laza)
        vals = vals1[0], vals1[1], vals2[0], vals2[1]

        legends = ["DAG/PI\n(rdga3)", "PA/PI\n(rdga3)", "DAG/PI\n(laza22)",
                   "PA/PI\n(laza22)"]
        ax.bar(ind, exp_observations, width, color='b', label="WT")
        ax.bar(ind + width, vals, width, color='r', label="MT")
        ax.set_xticks(ind + width / 2)
        plt.axhline(1, linestyle="--")
        plt.axhline(2.5, linestyle="--")
        ax.set_xticklabels(legends)
        plt.legend(loc=0)

    def set_layout(self):

        btn_layout = QVBoxLayout()
        btn_layout.addWidget(self.reset_btn)
        btn_layout.addWidget(self.print_btn)
        btn_layout.addWidget(self.toggle_plot)

        self.grid.addWidget(self.ss_label, 2, 1)
        self.grid.addWidget(self.toolbar, 1, 0)
        self.grid.addWidget(self.canvas, 2, 0)
        self.grid.addLayout(self.mutant_list, 3, 1)
        self.grid.addLayout(self.enzyme_parameters, 4, 0)
        self.grid.addLayout(btn_layout, 4, 1)
        self.setLayout(self.grid)


def interact():
    app = QApplication(sys.argv)
    main = MutantWindow()
    main.show()
    sys.exit(app.exec_())
