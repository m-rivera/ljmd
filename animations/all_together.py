from manim import *
import numpy as np
from numpy.random import default_rng

# colors
teal = "#3E969A"
l_teal = "#64e8ef"
pink = "#F24680"
yellow = "#FFA34F"
l_yellow = "#fcbc80"
green = "#3C9D00"
brown = "#87402f"
black = "#000000"
bg_col = black
grey = "#8e8e8e"

# params
sigma = 2
epsilon = 1
x_min_ener = 2 ** (1 / 6) * sigma


def cool_dot(pos):
    return Dot(np.array([pos[0],pos[1],0]), radius=0.15, color=teal, fill_opacity=1.)


def cool_arrow(start, end, col):
    return Arrow(
        np.array([start[0],start[1],0]), np.array([end[0],end[1],0]), color=col, stroke_width=5, buff=0
    )

def read_xyz(in_name):
    """
    Read a .xyz file.

    Works for files containing several configurations e.g. a relaxation
    trajectory.

    """
    with open(in_name) as xyz_file:
        xyz_content = xyz_file.readlines()

    # main list where each element is a relaxation step
    atom_step = []

    for i, line in enumerate(xyz_content):

        # if the line is the amount of atoms in the system
        if line.strip():
            if line.split()[0].isdigit():
                positions = []

                # from 2 lines after the amount of atoms to the last atom line
                # for the relaxation step
                for line_in_step in xyz_content[i + 2:i + int(line) + 2]:
                    elemAtom = line_in_step.split()[0]
                    xAtom, yAtom, zAtom = map(float, line_in_step.split()[1:])
                    positions.append([xAtom/3, yAtom/3, zAtom/3])

                atom_step.append(np.array(positions))

    xyz_file.close()
    return atom_step



class AllTogether(MovingCameraScene):
    def construct(self):

        self.camera.background_color = bg_col
        # set up plot
        numberplane = NumberPlane(background_line_style = {"stroke_color":grey})

        self.camera.frame.save_state()
        self.play(self.camera.frame.animate.scale(1))
        supercell_1 = Square(side_length=4, color = l_yellow, fill_color = grey, fill_opacity = 0.2)
        supercell_1.move_to(np.array([-3.0,0.0,0.0]))
        supercell_2 = Square(side_length=4, color = l_yellow, fill_color = grey, fill_opacity = 0.2)
        supercell_2.move_to(np.array([3.0,0.0,0.0]))

        rectangle_text_1 = Rectangle().next_to(supercell_1,UP)
        rectangle_text_2 = Rectangle().next_to(supercell_2,UP)
        rectangle_text_step = Rectangle().move_to(3*UP)
        text_1 = Tex("Low T").next_to(supercell_1,UP)
        text_2 = Tex("High T").next_to(supercell_2,UP)
        step_text = Text("Step").move_to(UP*3)
        step_int = Integer(1).move_to(UP*2.2)



        traj_1 = read_xyz("trajectory_solid.xyz")
        traj_2 = read_xyz("trajectory_liquid.xyz")

        x_vt_1 = [ValueTracker(traj_1[0][i][0] - 5) for i in range(32)]
        y_vt_1 = [ValueTracker(traj_1[0][i][1] - 2) for i in range(32)]


        x_vt_2 = [ValueTracker(traj_2[0][i][0] + 1) for i in range(32)]
        y_vt_2 = [ValueTracker(traj_2[0][i][1] - 2) for i in range(32)]


        dots_1 = [always_redraw(lambda i=i: cool_dot(np.array([x_vt_1[i].get_value(),y_vt_1[i].get_value()]))) for i in range(32)]
        dots_2 = [always_redraw(lambda i=i: cool_dot(np.array([x_vt_2[i].get_value(),y_vt_2[i].get_value()]))) for i in range(32)]

#        self.play(Create(numberplane),
#                Create(rectangle_text_1),
#                Create(rectangle_text_2),
#                Create(rectangle_text_step),
#                Create(supercell_1),
#                *[Create(i) for i in dots_1],
#                Create(supercell_2),
#                *[Create(i) for i in dots_2],
#                Create(text_1),
#                Create(text_2),
#                Create(step_text),
#                Create(step_int))
#        self.wait()

        self.play(Create(supercell_1),
                *[Create(i) for i in dots_1],
                Create(supercell_2),
                *[Create(i) for i in dots_2],
                Create(text_1),
                Create(text_2),
                Create(step_text),
                Create(step_int))
        self.wait()


        for i in range(1,100):
            self.play(*[x_vt_1[j].animate.set_value(traj_1[i][j][0] - 5) for j in range(32)],
                    *[y_vt_1[j].animate.set_value(traj_1[i][j][1] - 2) for j in range(32)],
                    *[x_vt_2[j].animate.set_value(traj_2[i][j][0] + 1) for j in range(32)],
                    *[y_vt_2[j].animate.set_value(traj_2[i][j][1] - 2) for j in range(32)],
                    #step_int.animate.set_value(i+2).move_to(UP*2.2),
                    step_int.animate.become(Integer(i+1)).move_to(UP*2.2),
                    run_time=0.1,
                    rate_func=rate_functions.linear,
                    )
            step_int.become(Integer(i+1)).move_to(UP*2.2)
        self.wait()
        self.play(*[FadeOut(i) for i in dots_1],
                FadeOut(supercell_1),*[FadeOut(i) for i in dots_2],FadeOut(supercell_2),
                FadeOut(text_1),
                FadeOut(text_2),
                FadeOut(step_text),
                FadeOut(step_int))
        self.wait()




