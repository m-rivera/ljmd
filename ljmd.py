# /usr/bin/env python
"""
Simple implementation of a molecular dynamics code for a Lennard-Jones system

Inspiration is taken from a Fortran code for liquid Argon written by Alex
Pacheco on 30/01/2014 and modified by Devis Di Tommaso on 08/02/2017.  The atoms
start at FCC positions, are assigned a random velocity of magnitudes fulfilling
an initial temperature requirement, and then evolve in a cubic periodic simulation
cell according to the velocity Verlet algorithm. A scaling thermostat is applied at each step.

The units are Angstrom, seconds and kJ/mol. Note that the kJ/mol energy is that of the total cell, not per atom.
03/03/2021 - Miguel Rivera

"""
import numpy as np
import time
from itertools import combinations

### Functions

def pos_to_xyz(pos,file_name):
    with open(file_name, "a") as file_out:
        file_out.write(str(len(pos)) + "\n\n")
        # for line in pos*1e10:
        for line in pos:
            file_out.write("Ar {:9.3f} {:9.3f} {:9.3f}\n".format(line[0], line[1], line[2]))
    return


def interatomic_potential(inter_dist):
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
    epsilon = 0.99363  # depth of the well in kJ/mol
    sigma = 3.405  # LJ interaction = 0 distance in A
    ener = 4 * epsilon * ((sigma / inter_dist) ** 12 - (sigma / inter_dist) ** 6)

    return ener


def interatomic_force(inter_dist):
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
    epsilon = 0.99363  # depth of the well in kJ/mol
    sigma = 3.405  # LJ interaction = 0 distance in A
    force = (
        48
        * epsilon
        * ((sigma ** 12) / (inter_dist ** 13) - 0.5 * (sigma ** 6) / (inter_dist ** 7))
    )
    return force


def cancelled_linear_momentum(vel):
    """
    Set the linear momentum of a velocity array to 0

    Parameters
    ----------
    vel : N x 3 numpy array
        Initial velocities
    Returns
    -------
    vel_cancelled : N x 3 numpy array
        Velocities after cancellation

    """
    # velocity of the centre of mass
    com_vel = np.sum(vel, axis=0)
    # remove COM velocity/number of atoms from all atoms
    vel_cancelled = vel - com_vel / len(vel)

    return vel_cancelled


def temp_from_vel(vel):
    """Return temperature in T from a velocity array in A/s"""
    # Angstrom**2 to meter**2 conversion
    prefactor = 1e-20
    temperature = 2*(
        prefactor
        * np.sum(np.linalg.norm(vel, axis=1) ** 2)
        * atomic_mass
        / (3 * g_const * n_atoms)
    )
    return temperature


def scaled_velocities(vel, temperature):
    """
    Scale a velocity vector according to a certain temperature

    Parameters
    ----------
    vel : N x 3 numpy array
        Initial velocities
    temperature : float
        Temperature in K
    Returns
    -------
    vel_scaled : N x 3 numpy array
        Scaled velocities

    """
    temp_initial = temp_from_vel(vel)
    vel_scaled = vel * np.sqrt(temperature / temp_initial)
    return vel_scaled

def rdf(pos,dr):
    """
    Return a radial distribution

    Arguments
    ---------
    pos : N x 3 numpy array
        Cartesian coordinates
    dr : float
        Radius step size between 0 and box_side length in Angstrom
    Returns
    -------
    x_range : numpy array
        X coordinate of the beginning of each bin
    left_bin_edges : numpy array
        Normalised count of interatomic distances within the bin. The normalisation is by the volume of the shell.

    """
    x_range = np.arange(0, box_side+dr, dr)

    # all distances with no double counting
    comb = np.array(list(combinations(pos, 2)))
    dist_vecs = comb[:, 1, :] - comb[:, 0, :]
    dist_vecs_pbc = (dist_vecs + box_side / 2) % box_side - box_side / 2
    distances = np.linalg.norm(dist_vecs_pbc,axis=1)

    # sort counts in a histogram
    counts, tmp_range = np.histogram(distances, x_range)
    left_bin_edges = x_range[:-1]

    # normalise by number of distances and particle density
    norm_counts = 2*counts/(((left_bin_edges+dr)**3 - left_bin_edges**3)* 4*np.pi/3 * n_atoms**2/box_side**3)

    # only return first half of the data, the latter is ruined by PBC
    return left_bin_edges[:round(len(left_bin_edges)/2)], norm_counts[:round(len(norm_counts)/2)]

### Main code

# timer
start_time = time.perf_counter()

# parameters
n_cell = 4 # number of unit cells along one axis of the supercell
lattice_const = 5.256  # in A
box_side = n_cell * lattice_const  # side of the simulation box in A
n_steps = 50 # number of simulation steps
temp = 87.81  # in K
time_step = 2e-11  # s
atomic_mass = 39.95  # in u
g_const = 8.31446261815324e-3  # gas constant in kJ / (mol * K)
pos_file = "trajectory.xyz"
report_file = "report.csv"

### Initial positions

# cubic lattice
lattice_vectors = np.array(
    [
        [lattice_const, 0.0, 0.0],
        [0.0, lattice_const, 0.0],
        [0.0, 0.0, lattice_const],
    ]
)

# one FCC cell
fcc_pos = (
    np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
    * lattice_const
)

# supercell of FCC cells
pos_i = []
for i in range(n_cell):
    for j in range(n_cell):
        for k in range(n_cell):
            # generate a translation vector
            trans_mult = np.array([i, j, k])
            trans_vec = np.sum(lattice_vectors * trans_mult, axis=1)

            # append new cell
            new_cell = fcc_pos + trans_vec
            pos_i.extend(new_cell.tolist())

pos_i = np.array(pos_i)

### Assign starting velocities

# gaussian distribution number generator
n_atoms = (
    len(fcc_pos) * n_cell ** 3
)

# random seed
rng = np.random.default_rng(123456789)
vel_i = rng.standard_normal((n_atoms, 3))

# set linear momentum to 0
vel_i = cancelled_linear_momentum(vel_i)

# scale velocities according to temperature
vel_i = scaled_velocities(vel_i, temp)

### Simulation loop

# start with null acceleration
acc_i = np.zeros((n_atoms, 3))

# start a output files
with open(pos_file, "w") as fname:
    pos_to_xyz(pos_i,pos_file)

with open(report_file,"w") as fname:
    fname.write("Step, Temperature, Energy\n")

# loop each MD step
for i in range(n_steps):
    # velocity Verlet algorithm position update
    pos_f = pos_i + vel_i * time_step + (acc_i * time_step ** 2) / 2

    # apply periodic boundary conditions
    pos_f = pos_f % box_side

    # interatomic distances
    ener = 0
    force = np.zeros((n_atoms, 3))

    ## Apply potentials to the positions using numpy vectorisation
    # all pair positions
    comb = np.array(list(combinations(pos_f, 2)))

    # distance vector between each pair
    dist_vecs = comb[:, 1, :] - comb[:, 0, :]

    # periodic boundary
    dist_vecs_pbc = (dist_vecs + box_side / 2) % box_side - box_side / 2

    # array of scalar distances
    distances = np.linalg.norm(dist_vecs_pbc,axis=1)

    # calculate energies
    ener = np.sum(interatomic_potential(distances))

    # array of vector force per pair
    force_by_pair = dist_vecs_pbc / distances[:,None] * np.array(interatomic_force(distances))[:,None]

    # reshape it into a n_atoms x n_atoms (x 3) upper matrix
    force_mat = np.zeros((n_atoms,n_atoms,3))
    mask = np.triu_indices(n_atoms,k=1)
    force_mat[mask] = force_by_pair

    # add up rows and take away columns for a total
    force = np.sum(force_mat,axis=0) - np.sum(force_mat,axis=1)

    ## New velocities form forces
    # F = ma with 1e20 for Angstrom**2 to meter**2
    acc_f = force * 1e20 / atomic_mass

    # velocity update
    vel_f = vel_i + (acc_f + acc_i) / 2 * time_step

    # set linear momentum to 0
    vel_f = cancelled_linear_momentum(vel_f)

    # scale velocities according to temperature
    vel_f = scaled_velocities(vel_f, temp)

    # set up for next step
    pos_i = pos_f
    vel_i = vel_f
    acc_i = acc_f

    # print out report
    print(
        "Step : {:>4} Temperature: {:>7.3f} K Cell energy: {:>9.3f} kJ/mol".format(
            i + 1, temp_from_vel(vel_f), ener
        )
    )

    with open(report_file,"a") as fname:
        fname.write("{},{},{}\n".format(i+1,temp_from_vel(vel_f),ener))

    # print out new goemetry
    pos_to_xyz(pos_f, "trajectory.xyz")

# timer
end_time = time.perf_counter()
print("Run time: {:>8.4f} s.".format(end_time - start_time))

rdf_distance, radial_density = rdf(pos_f,.1)

import matplotlib.pyplot as plt
plt.plot(rdf_distance, radial_density)
plt.savefig("rdf.png")

