from manim import *
import numpy as np

# colors
teal = "#3E969A"
pink = "#F24680"
yellow = "#FFA34F"
green = "#3C9D00"
brown = "#87402f"
bg_col = "#000000"
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

class Verlet(Scene):
    def construct(self):
        self.camera.background_color = bg_col
        # number line
        num_line = NumberLine()
        num_line.shift(DOWN * 2.5)
        self.play(Create(num_line))
        num_line_label = Tex(r"r")
        num_line_label.shift(DOWN * 3)
        self.play(Write(num_line_label))

        # atoms
        r_vt = ValueTracker(-4)
        atom_i = always_redraw(
            lambda: cool_dot(num_line.number_to_point(r_vt.get_value()) + UP)
        )
        atom_i_lab = always_redraw(lambda: MathTex(r"i").move_to(atom_i))
        atom_j = always_redraw(
            lambda: cool_dot(num_line.number_to_point(-1 * r_vt.get_value()) + UP)
        )
        atom_j_lab = always_redraw(lambda: MathTex(r"j").move_to(atom_j))

        # r label
        r_lab = MathTex(r"\mathbf{r}^i_n")
        r_lab.set_color(teal)
        r_lab.shift(UP * 3 + LEFT)
        self.play(
            Create(atom_i),
            Write(atom_i_lab),
            Create(atom_j),
            Write(atom_j_lab),
            Write(r_lab),
        )
        self.wait(1)

        # v label
        v_vt = ValueTracker(1.8)
        v_lab = MathTex(r"\mathbf{v}^i_n")
        v_lab.set_color(yellow)
        v_lab.shift(2.2 * UP + LEFT)
        v_i = always_redraw(
            lambda: cool_arrow_pos(
                v_vt.get_value(),
                num_line.number_to_point(r_vt.get_value()) + 1.8 * UP,
                yellow,
            )
        )
        v_j = always_redraw(
            lambda: cool_arrow_neg(
                v_vt.get_value(),
                num_line.number_to_point(-1 * r_vt.get_value()) + 1.8 * UP,
                yellow,
            )
        )
        self.play(Create(v_i), Create(v_j), Write(v_lab))

        self.wait(1)

        # n+1 labels
        arrow_comp = MathTex(r"\rightarrow")
        arrow_comp.shift(2.5 * UP)
        delta_t_lab = MathTex(r"\delta t")
        delta_t_lab.shift(3.0 * UP)
        r_lab_next = MathTex(r"\mathbf{r}^i_{n+1}")
        v_lab_next = MathTex(r"\mathbf{v}^i_{n+1}")
        r_lab_next.set_color(teal)
        v_lab_next.set_color(yellow)
        r_lab_next.shift(UP * 3 + RIGHT)
        v_lab_next.shift(UP * 2.2 + RIGHT)
        self.play(
            Write(arrow_comp), Write(delta_t_lab), Write(v_lab_next), Write(r_lab_next)
        )
        self.play(v_vt.animate.set_value(2.0), r_vt.animate.set_value(-3))


        self.play(v_vt.animate.set_value(1.8), r_vt.animate.set_value(-4))

        self.wait(1)

        self.play(
            FadeOut(r_lab),
            FadeOut(v_lab),
            FadeOut(r_lab_next),
            FadeOut(v_lab_next),
            FadeOut(arrow_comp),
            FadeOut(delta_t_lab),
        )

        # Verlet equations

        # r Verlet
        r_eq = MathTex(
            r"\mathbf{r}^i_{n+1}",
            r"= ",
            r"\mathbf{r}^i_n",
            r" + ",
            r"\mathbf{v}^i_n \delta t",
            r"+ ",
            r"\frac{1}{2} \mathbf{a}^i_n \delta t^2",
        )
        r_eq[0].set_color(teal)
        r_eq[2].set_color(teal)
        r_eq[4][0:3].set_color(yellow)
        r_eq[6][3:6].set_color(pink)
        r_eq.shift(UP * 3)

        self.play(Write(r_eq))
        self.wait(1)


        # v Verlet
        v_eq = MathTex(
            r"\mathbf{v}^i_{n+1}",
            r" = ",
            r"\mathbf{v}^i_n",
            r" + ",
            r"\frac{\mathbf{a}^i_n + \mathbf{a}^i_{n+1}}{2} \delta t",
        )
        v_eq[0].set_color(yellow)
        v_eq[2].set_color(yellow)
        v_eq[4][0:3].set_color(pink)
        v_eq[4][4:9].set_color(pink)
        v_eq.shift(UP * 1.8)

        self.play(Write(v_eq))

        self.wait()

        # force
        self.play(r_eq.animate.scale(0.8), v_eq.animate.scale(0.8))
        a_eq = MathTex(
            r"\mathbf{a}^i_{n+1}",
            r" =",
            r"\frac{\sum_{j \neq i} \mathbf{F}_{LJ}(\abs{\mathbf{r}^i_{n+1}-\mathbf{r}^j_{n+1})}}{m_i}",
        )
        a_eq[0].set_color(pink)
        a_eq[2][12:17].set_color(teal)
        a_eq.move_to(v_eq)
        a_eq.scale(0.7)
        self.play(v_eq.animate.shift(1.2 * DOWN))
        self.wait()
        a_vt = ValueTracker(0.5)
        a_i = always_redraw(
            lambda: cool_arrow_pos(
                a_vt.get_value(),
                num_line.number_to_point(r_vt.get_value()) + 2 * UP,
                pink,
            )
        )
        a_j = always_redraw(
            lambda: cool_arrow_neg(
                a_vt.get_value(),
                num_line.number_to_point(-1 * r_vt.get_value()) + 2 * UP,
                pink,
            )
        )
        self.play(Write(a_eq), Create(a_i), Create(a_j))
        self.wait()

        # algorithm
        shift = RIGHT * 3.5
        self.play(
            r_eq.animate.shift(shift),
            v_eq.animate.shift(shift),
            a_eq.animate.shift(shift),
        )
        scale = 0.8
        n_vt = ValueTracker(0)
        shift = UP * 3.2 + LEFT * 5.5
        algo_text = [
            always_redraw(
                lambda: Tex(r"$n=" + str(int(n_vt.get_value())) + r"$").shift(shift)
            ),
            Tex(r"1. Calculate the positions at $n+1$.").scale(scale),
            Tex(r"2. Calculate the acceleration at $n+1$.").scale(scale),
            Tex(r"3. Calculate the velocity at $n+1$.").scale(scale),
            Tex(r"4. Increase $n$ by 1.").scale(scale),
            Tex(r"5. Loop back to step 1.").scale(scale),
        ]
        # shift = UP
        for i, text in enumerate(algo_text[1:]):
            text.move_to(algo_text[0])
            text.align_to(algo_text[0], LEFT)
            text.shift((i+1)*DOWN / 2)

        self.play(*[Write(i) for i in algo_text])
        self.wait()

        # step 1
        boxes_algo = [
            SurroundingRectangle(i, buff=0.1, color=grey) for i in algo_text[1:]
        ]
        box_r_i = always_redraw(
            lambda: SurroundingRectangle(atom_i, buff=0.1, color=grey)
        )
        box_r_eq = SurroundingRectangle(r_eq, buff=0.1, color=grey)
        self.play(
            Create(boxes_algo[0]), Create(box_r_i), Create(box_r_eq)
        )
        self.wait()
        self.play(r_vt.animate.set_value(-2.5))
        self.wait()

        # step 2
        box_a_i = always_redraw(lambda: SurroundingRectangle(a_i, buff=0.1, color=grey))
        box_a_eq = SurroundingRectangle(a_eq, buff=0.1, color=grey)
        self.play(
            ReplacementTransform(boxes_algo[0], boxes_algo[1]),
            ReplacementTransform(box_r_i, box_a_i),
            ReplacementTransform(box_r_eq, box_a_eq),
        )
        self.wait()
        self.play(a_vt.animate.set_value(0.7))
        self.wait()

        # step 3
        box_v_i = always_redraw(lambda: SurroundingRectangle(v_i, buff=0.1, color=grey))
        box_v_eq = SurroundingRectangle(v_eq, buff=0.1, color=grey)
        self.play(
            ReplacementTransform(boxes_algo[1], boxes_algo[2]),
            ReplacementTransform(box_a_i, box_v_i),
            ReplacementTransform(box_a_eq, box_v_eq),
        )
        self.wait()
        self.play(v_vt.animate.set_value(2.0))
        self.wait()

        # step 4
        box_n = SurroundingRectangle(algo_text[0], buff=0.1, color=grey)
        self.play(
            ReplacementTransform(boxes_algo[2], boxes_algo[3]),
            ReplacementTransform(box_v_i, box_n),
            FadeOut(box_v_eq),
        )
        self.wait()
        self.play(n_vt.animate.set_value(1))
        self.wait()

        # step 5
        self.play(FadeOut(box_n), ReplacementTransform(boxes_algo[3], boxes_algo[4]))
        self.wait()
        last_box = SurroundingRectangle(algo_text[1], buff=0.1, color=grey)
        self.play(ReplacementTransform(boxes_algo[4], last_box))
        self.wait()

        # move
        self.play(
            FadeOut(last_box),
            FadeOut(r_eq),
            FadeOut(a_eq),
            FadeOut(v_eq),
            FadeOut(algo_text[1]),
            FadeOut(algo_text[2]),
            FadeOut(algo_text[3]),
            FadeOut(algo_text[4]),
            FadeOut(algo_text[5]),
            FadeOut(atom_i),
            FadeOut(atom_j),
            FadeOut(v_i),
            FadeOut(v_j),
            FadeOut(a_i),
            FadeOut(a_j),
            FadeOut(num_line),
            FadeOut(num_line_label),
            FadeOut(atom_i_lab),
            FadeOut(atom_j_lab),
            FadeOut(algo_text[0]),
        )

        self.wait()
