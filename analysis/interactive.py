import sys
from functools import partial

import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QApplication, QGridLayout, \
    QLabel, QVBoxLayout, QRadioButton, QSlider, QButtonGroup, QWidget, \
    QComboBox, QHBoxLayout, QCheckBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import \
    NavigationToolbar2QT as NavigationToolbar

from analysis.helper import *
from analysis.ode_analysis import ODEAnalysis

LIGHT_STIMULATION = 85
CARRY_FACTORS = ["x 1", "x 10", "x 0.1", "x 0.01"]
PARAMETER_NO = 0


def provide_condition(base):
    base.buffer_time = 0.5
    base.normalize_enzymes()
    base.initialize()
    base.stimulate(LIGHT_STIMULATION)
    base.recover(2)


class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.setWindowTitle('Feedback Interactive Analysis')

        self.enzyme = E_PIP5K
        self.f_type = FEEDBACK_POSITIVE
        self.f_carry = 1
        self.f_multi = 1
        self.f_hill = 1
        self.f_sub = I_PIP2
        self.carry_factor = 1
        self.lipids_to_analyze = [I_PIP2]

        # a figure instance to plot on
        self.figure = plt.figure()
        # this is the Canvas Widget that displays the `figure`
        self.canvas = FigureCanvas(self.figure)
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        grid = QGridLayout()
        widget = QWidget(self)  # central widget

        s_box = QHBoxLayout()

        substrate_group = QButtonGroup(widget)  # Substrate group
        s_box.addStretch()
        s_box.setContentsMargins(1, 15, 1, 1)
        s_box.addWidget(QLabel("Feedback From"))
        self.substrates = []
        for l in LIPID_NAMES:
            b = QRadioButton(l, self)
            b.setObjectName(l)
            b.setCheckable(True)
            if l == L_PIP2:
                b.setChecked(True)
            b.clicked.connect(partial(self.substrate_clicked, l))
            self.substrates.append(b)
            substrate_group.addButton(b)
            s_box.addWidget(b)

        display_box = QVBoxLayout()
        display_box.setContentsMargins(30, 1, 1, 1)
        display_box.addStretch()
        self.display_lipids = []
        for l in LIPID_NAMES:
            b = QCheckBox(l)
            b.setObjectName(l)
            b.setCheckable(True)
            if l == L_PIP2:
                b.setChecked(True)
            b.clicked.connect(partial(self.on_display_lipid_changed, l))
            self.display_lipids.append(b)
            display_box.addWidget(b)
        display_box.addStretch()

        feedback_group = QButtonGroup(widget)  # Number group
        f_box = QVBoxLayout()
        f_box.setContentsMargins(30, 1, 1, 1)
        f_box.addStretch()
        f_box.addWidget(QLabel("Feedback Type"))
        self.b_pos_feed = QRadioButton("Positive", self)
        self.b_pos_feed.setObjectName("p")
        self.b_pos_feed.setCheckable(True)
        self.b_pos_feed.setChecked(True)
        self.b_pos_feed.clicked.connect(partial(self.feedback_changed, "p"))
        feedback_group.addButton(self.b_pos_feed)
        f_box.addWidget(self.b_pos_feed)

        self.b_neg_feed = QRadioButton("Negative", self)
        self.b_neg_feed.setObjectName("n")
        self.b_neg_feed.setCheckable(True)
        self.b_neg_feed.clicked.connect(partial(self.feedback_changed, "n"))
        feedback_group.addButton(self.b_neg_feed)
        f_box.addWidget(self.b_neg_feed)
        f_box.addStretch()

        a_box = QVBoxLayout()
        a_box.setContentsMargins(30, 30, 30, 15)
        a_slider = QSlider(Qt.Horizontal, self)
        a_slider.setMinimum(1)
        a_slider.setMaximum(100)
        a_slider.setTickInterval(1)
        self.a_label = QLabel("Multiplication Factor : 1")
        a_slider.valueChanged[int].connect(self.multi_changed)
        a_slider.setValue(1)
        a_box.addWidget(self.a_label)
        a_box.addWidget(a_slider)

        c_box = QVBoxLayout()
        c_box.setContentsMargins(30, 0, 30, 15)
        c_slider = QSlider(Qt.Horizontal, self)
        c_slider.setMinimum(1)
        c_slider.setMaximum(10)
        c_slider.setTickInterval(1)
        self.c_label = QLabel("Carrying Capacity : 1")
        c_slider.valueChanged[int].connect(self.carry_changed)
        c_slider.setValue(1)

        c_sub_box = QHBoxLayout()
        c_sub_box.addWidget(self.c_label)

        carry_factor = QComboBox(self)
        for c in CARRY_FACTORS:
            carry_factor.addItem(c)
        carry_factor.activated[str].connect(self.on_carry_factor)
        c_sub_box.addWidget(carry_factor)
        c_sub_box.addStretch()

        c_box.addLayout(c_sub_box, 1)
        c_box.addWidget(c_slider)

        h_box = QVBoxLayout()
        h_box.setContentsMargins(30, 0, 30, 15)
        h_slider = QSlider(Qt.Horizontal, self)
        h_slider.setMinimum(0)
        h_slider.setMaximum(5)
        h_slider.setTickInterval(1)
        self.h_label = QLabel("Hill Coefficient")
        h_slider.valueChanged[int].connect(self.hill_changed)
        h_slider.setValue(1)
        h_box.addWidget(self.h_label)
        h_box.addWidget(h_slider)

        e_box = QVBoxLayout()
        e_box.setContentsMargins(30, 1, 1, 1)
        e_box.addStretch()
        e_box.addWidget(QLabel("Feedback to"))
        enz_box = QComboBox(self)
        for e in ENZYME_NAMES:
            enz_box.addItem(e)
        enz_box.activated[str].connect(self.on_enzyme_selection)
        e_box.addWidget(enz_box)
        e_box.addStretch()

        grid.addWidget(self.toolbar, 1, 0)
        grid.addWidget(self.canvas, 2, 0)
        grid.addLayout(display_box, 1, 1, 2, 2)
        grid.addLayout(s_box, 4, 0)
        grid.addLayout(f_box, 5, 1, 2, 2)
        grid.addLayout(a_box, 5, 0)
        grid.addLayout(c_box, 6, 0)
        grid.addLayout(h_box, 7, 0)
        grid.addLayout(e_box, 7, 1)
        self.setLayout(grid)

    def plot(self):

        # instead of ax.hold(False)
        self.figure.clear()

        enz = extract_enzyme_set("para.txt", PARAMETER_NO)
        enz2 = extract_enzyme_set("para.txt",
                                  PARAMETER_NO)  # Required because ODEAnalysis
        # changes enzyme values
        base = ODEAnalysis(S_OPEN_2, 1, enz)
        feed = ODEAnalysis(S_OPEN_2, 1, enz2, self.feed_para)

        provide_condition(base)
        provide_condition(feed)

        ax = self.figure.add_subplot(111)

        for lipid in self.lipids_to_analyze:
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

        # refresh canvas
        self.canvas.draw()

    @property
    def feed_para(self):
        return {
            self.enzyme: {
                F_TYPE_OF_FEEDBACK: self.f_type,
                F_HILL_COEFFICIENT: self.f_hill,
                F_MULTIPLICATION_FACTOR: self.f_multi,
                F_CARRYING_CAPACITY: self.f_carry * self.carry_factor,
                F_FEED_SUBSTRATE_INDEX: self.f_sub
            }
        }

    def on_carry_factor(self, text):
        if text == CARRY_FACTORS[0]:
            self.carry_factor = 1
        elif text == CARRY_FACTORS[1]:
            self.carry_factor = 10
        elif text == CARRY_FACTORS[2]:
            self.carry_factor = 0.1
        elif text == CARRY_FACTORS[3]:
            self.carry_factor = 0.01

        self.plot()

    def on_enzyme_selection(self, text):
        self.enzyme = text
        self.plot()

    def on_display_lipid_changed(self, text):
        self.lipids_to_analyze.clear()
        for b in self.display_lipids:
            if b.isChecked():
                self.lipids_to_analyze.append(LIPID_NAMES.index(b.objectName()))

        self.plot()

    def substrate_clicked(self, e):
        for b in self.substrates:
            bt = b  # type: QRadioButton
            if bt.objectName() != e:
                if bt.isChecked():
                    bt.toggle()
        self.f_sub = LIPID_NAMES.index(e)
        self.plot()

    def multi_changed(self, value):
        multi_value = value
        multi_label = "Multiplication Factor : %.2f" % multi_value
        self.a_label.setText(multi_label)
        self.f_multi = value
        self.plot()

    def carry_changed(self, value):
        carry_value = value
        carry_label = "Carrying capacity : %.2f" % carry_value
        self.c_label.setText(carry_label)
        self.f_carry = value
        self.plot()

    def hill_changed(self, value):
        hill_value = value
        hill_label = "Hill Coefficient : %.2f" % hill_value
        self.h_label.setText(hill_label)
        self.f_hill = value
        self.plot()

    def feedback_changed(self, name):
        if name == "p":
            self.f_type = FEEDBACK_POSITIVE
            self.b_neg_feed.setChecked(False)
        else:
            self.f_type = FEEDBACK_NEGATIVE
            self.b_neg_feed.setChecked(False)

        self.plot()


def interact():
    app = QApplication(sys.argv)
    main = Window()
    main.show()
    sys.exit(app.exec_())
