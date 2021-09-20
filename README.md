# Lennard-Jones Molecular Dynamics

## Details

This is an educational Python code for Lennard-Jones Molecular Dynamics. It was
adapted from a similar Fortran code written by Alex Pacheco on 30/01/2014, and
modified by Devis Di Tommaso on 08/02/2017.

It implements:

- The velocity Verlet algorithm.
- FCC initial conditions.
- A velocity scaling thermostat.
- Periodic boundary conditions with the minimal image convention.
- Radial distribution function generator.

There is an additional Jupyter notebook which contains the same code with blank
spaces for students to annotate.
The `animations/` contains files to make instructive animations using `manim`
v10.0.0. They are beginner animations, so don't expect best practices.

## Outputs

The program generates out four files:

- `trajectory.xyz`: The positions of the atoms.
- `report.csv`: The temperature and energy at each step.
- `rdf.png`: A rough figure of a radial distribution function (no labels).

The outputs are deliberately terse in order to encourage the learner to build
upon them.

## Dependencies

- numpy: For vector arithmetic.
- matplotlib: For plotting

## Comments

The code is layed out in a roughly efficient manner, but the objective is to be
pedagogical and have it easily modifiable by a learner. The only exception to
this is lines 256 to 281, which uses numpy vectorisation (not a beginner
friendly technique).

An explicit for-loop version would be easier to write and understand, but the
increase in efficiency with numpy is an invaluable improvement in the quality of
life of the user.
