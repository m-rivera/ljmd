from manim import *
import numpy as np

# colors
teal = "#3E969A"
pink = "#F24680"
yellow = "#FFA34F"
green = "#3C9D00"
l_green = "#87402f"
l_green = "#E2F89C"
d_blue = "#032B43"
purple = "#55286F"
d_purple = "#210B2C"
black = "#000000"
bg_col = black
grey = "#8e8e8e"

# params
sigma = 2
epsilon = 1
x_min_ener = 2 ** (1 / 6) * sigma


def cool_dot(pos):
    return Dot(pos, radius=0.5, color=teal, fill_opacity=0.5)


def cool_arrow_pos(mag, origin, col):
    return Arrow(
        origin, origin + np.array([mag, 0, 0]), color=col, stroke_width=5, buff=0
    )


def cool_arrow_neg(mag, origin, col):
    return Arrow(
        origin, origin - np.array([mag, 0, 0]), color=col, stroke_width=5, buff=0
    )


def make_pot_func(sigma, epsilon):
    return (
        lambda inter_dist: 4
        * epsilon
        * ((sigma / inter_dist) ** 12 - (sigma / inter_dist) ** 6)
    )


def interatomic_potential_generic(inter_dist, sigma=1.2, epsilon=1):
    """
    Return the energy of the van der Waals interaction given an interatomic distance

    Parameters
    ----------
    inter_dist : float
        Distance between two atoms in Angstroms
    Returns
    -------
    ener : float
        Interaction energy in kJ / mol

    """
    ener = 4 * epsilon * ((sigma / inter_dist) ** 12 - (sigma / inter_dist) ** 6)

    return ener


def make_for_func(sigma, epsilon):
    return lambda inter_dist: (
        48
        * epsilon
        * ((sigma ** 12) / (inter_dist ** 13) - 0.5 * (sigma ** 6) / (inter_dist ** 7))
    )


def interatomic_force_generic(inter_dist, sigma=1.2, epsilon=1):
    """
    Return the force resulting from the interatomic van der Waals interaction

    Parameters
    ----------
    inter_dist : float
        Distance between two atoms in Angstroms
    Returns
    -------
    force : float
        Magnitude of the resulting force in kJ / (mol * Angstrom)

    """
    force = (
        48
        * epsilon
        * ((sigma ** 12) / (inter_dist ** 13) - 0.5 * (sigma ** 6) / (inter_dist ** 7))
    )
    return force


def interatomic_potential(inter_dist):
    return interatomic_potential_generic(inter_dist, sigma, epsilon)


def interatomic_potential_tight(inter_dist):
    return interatomic_potential_generic(inter_dist, sigma_tight, epsilon_tight)


def interatomic_potential_wide(inter_dist):
    return interatomic_potential_generic(inter_dist, sigma_wide, epsilon_wide)


class LJPot(Scene):
   # def __init__(self, **kwargs):
   #     GraphScene.__init__(
   #         self,
   #         x_min=0,
   #         x_max=5,
   #         # num_graph_anchor_points=100,
   #         y_min=-1.5,
   #         y_max=1.5,
   #         y_axis_label=r"$V_{LJ}(r), \abs{F_{LJ}}$",
   #         x_axis_label=r"$r$",
   #         graph_origin=LEFT * 5,
   #         **kwargs
   #     )

    def construct(self):
        ax = Axes(x_range = (0,5),
                y_range = (-1.5,1.5),
                axis_config={"include_tip": False})

        labels = ax.get_axis_labels(x_label="r", y_label=r"V_{LJ}(r)")

        self.camera.background_color = bg_col
        # equation
        lj_eq = MathTex(
            r"V_{LJ}(r) = 4\epsilon[(\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^6]"
        )
        lj_eq[0][0:3].set_color(teal)
        lj_eq[0][8].set_color(yellow)
        lj_eq[0][11].set_color(pink)
        lj_eq[0][19].set_color(pink)
        self.play(Write(lj_eq))
        self.wait(1)
        label_eps = Tex(r"$\epsilon$: Well depth")
        label_eps[0][0].set_color(yellow)
        label_sig = Tex(r"$\sigma$: Particle size")
        label_sig[0][0].set_color(pink)
        label_eps.shift(DOWN)
        label_sig.shift(DOWN * 1.5)
        self.play(Write(label_eps))
        self.play(Write(label_sig))
        self.wait(1)

        shift = UP * 3 + RIGHT * 3
        self.play(
            ApplyMethod(lj_eq.shift, shift),
            ApplyMethod(label_eps.shift, shift),
            ApplyMethod(label_sig.shift, shift),
        )

        # graph ener
        x_left = 1.13
        x_right = 5
        x_min_ener = 2 ** (1 / 6) * sigma

        self.play(Create(ax),Create(labels))

        vt_sigma = ValueTracker(sigma)
        vt_epsilon = ValueTracker(epsilon)

        graph = always_redraw(
            lambda: ax.get_graph(
                make_pot_func(vt_sigma.get_value(), vt_epsilon.get_value()),
                x_range=[x_left,x_right],
                #x_min=x_left,
                #x_max=x_right,
                color=teal,
            )
        )
        graph.set_stroke(width=7, color=teal)
        self.add(graph)
        self.play(Create(graph), run_time=3)
        self.wait(1)

        # parameter lines
        eps_line = always_redraw(
            lambda: Line(
                ax.coords_to_point(2 ** (1 / 6) * vt_sigma.get_value(), 0),
                ax.coords_to_point(
                    2 ** (1 / 6) * vt_sigma.get_value(), -1 * vt_epsilon.get_value()
                ),
                color=yellow,
                #dash_length = 0.2,
                #dashed_ratio = 0.5,
                stroke_width=7,
                stroke_opacity=0.6
            )
        )
        sig_line = always_redraw(
            lambda: Line(
                ax.coords_to_point(0, 0),
                ax.coords_to_point(vt_sigma.get_value(), 0),
                color=pink,
                #dash_length = 0.2,
                #dashed_ratio = 0.5,
                stroke_width=7,
                stroke_opacity=0.6
            )
        )
        # bring back graph to front
        self.bring_to_front(graph)
        self.play(Create(eps_line), Create(sig_line))
        self.wait(1)

        # fade out params
        self.play(FadeOut(label_eps), FadeOut(label_sig))
        self.wait(1)

        # force
        lj_force = MathTex(
            r"|F_{LJ}(r)| = 48 \epsilon [ \frac{\sigma^{12}}{r^{13}} - \frac{\sigma^{6}}{2r^{7}}]"
        )
        lj_force[0][1:4].set_color(l_green)
        lj_force[0][11].set_color(yellow)
        lj_force[0][13].set_color(pink)
        lj_force[0][21].set_color(pink)
        lj_force.align_to(lj_eq)
        lj_force.shift(shift + DOWN)
        self.play(Write(lj_force))
        self.wait(1)

        graph_for = always_redraw(
            lambda: ax.get_graph(
                make_for_func(vt_sigma.get_value(), vt_epsilon.get_value()),
                x_range = [x_left,x_right],
                color=l_green,
            )
        )
        graph_for.set_stroke(width=7, color=l_green)
        self.add(graph_for)
        self.play(Create(graph_for), run_time=3)
        new_labels = ax.get_axis_labels(x_label="r", y_label=r"V_{LJ}(r), \abs{F_{LJ}}")
        self.play(Transform(labels,new_labels))
        self.wait(1)

        ## morph
        self.play(vt_sigma.animate.set_value(1.8), vt_epsilon.animate.set_value(1.2))
        self.wait(1)
        self.play(vt_sigma.animate.set_value(1.2), vt_epsilon.animate.set_value(1.2))
        self.wait(1)
        self.play(
            vt_sigma.animate.set_value(sigma), vt_epsilon.animate.set_value(epsilon)
        )
        self.wait(1)

        # FadeOut
        self.play(
            FadeOut(graph),
            FadeOut(graph_for),
            FadeOut(ax),
            FadeOut(lj_eq),
            FadeOut(lj_force),
            FadeOut(sig_line),
            FadeOut(eps_line),
            FadeOut(labels)
        )
        self.wait(1)

