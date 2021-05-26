#!/usr/bin/env python3

# you can install mendeleev package by typing 'pip3 install mendeleev' in terminal
# more informations https://pypi.org/project/mendeleev/#modal-close
from mendeleev import element
import numpy as np
import math
from itertools import islice
import sys
import os.path


# function for yes/no questions
def query_yes_no(question, default="no"):
    global prompt
    valid = {"yes": True, "y": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Wrong answer. Respond with 'yes' or 'no'.\n")


# vector normalisation
def unit_vector(vector):
    return vector / np.linalg.norm(vector)


# calculating angle between vectors in degrees
def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))


# cell parameters
class Cell:
    def __init__(self):
        self.a = []
        self.b = []
        self.c = []
        self.alpha = []
        self.beta = []
        self.gamma = []


freqs = []
create_movie = False

if len(sys.argv) < 2:
    print("You are in manual mode. \nTo load input data from file restart this program and specify the input file as: \n\"python phonons.py input_file.txt\"")

    # getting name of a vasp files
    file_name = input("To exit program press enter or specify the name (without extension) of the vasp files: ")
    if file_name == "":
        sys.exit(0)

    # frequency input
    while True:
        try:
            freq_min, freq_max = (float(s) for s in input('Enter the frequency range (in cm-1) you want to analyze '
                                                          '(with space as separator, e.g., 22 27): ').split())
        except ValueError:
            print('Incorrect data format. Try again: ')
        else:
            break
    
    # amplitude input
    while True:
        try:
            amp_min, amp_max, amp_step = (float(s) for s in input('Enter the amplitude range you want to analyze'
                                                                        ' with number of steps'
                                                                        ' (with space as separator,'
                                                                        ' e.g., -2 2 0.5): ').split())
        except ValueError:
            print('Incorrect data format. Try again: ')
        else:
            break
        amp_range = input("")

    # movies
    create_movie = query_yes_no('Do you want to create a .pdb file with phonon visualisation?')

else:
    with open(sys.argv[1], 'r') as input_file:
        for i, line in enumerate(input_file):
            if i == 1:
                file_name = str(line).rstrip()
            elif i == 3:
                amp_min = float(line.split()[0])
                amp_max = float(line.split()[1])
                amp_step = float(line.split()[2])
            elif i == 5:
                if str(line).rstrip() == 'True':
                    create_movie = True
            elif i == 7:
                line = line.rstrip().split()
                for freq in line:
                    freqs.append(str(freq))
                break

# checking if all files are in directory
poscar_file = file_name + '.poscar'
outcar_file = file_name + '.outcar'
if os.path.isfile(poscar_file) is False:
    sys.exit("Missing file " + poscar_file + ". Check it and try again.")
if os.path.isfile(outcar_file) is False:
    sys.exit("Missing file " + outcar_file + ". Check it and try again.")

cell_array = np.zeros((3, 3))

# read first POSCAR file lines
with open(poscar_file, 'r') as pos_file:
    for i, line in enumerate(pos_file):
        if i == 1:
            lattice_constant = float(line.split()[0])
        elif i == 2:
            cell_array[0, 0] = float(line.split()[0])
            cell_array[0, 1] = float(line.split()[1])
            cell_array[0, 2] = float(line.split()[2])
        elif i == 3:
            cell_array[1, 0] = float(line.split()[0])
            cell_array[1, 1] = float(line.split()[1])
            cell_array[1, 2] = float(line.split()[2])
        elif i == 4:
            cell_array[2, 0] = float(line.split()[0])
            cell_array[2, 1] = float(line.split()[1])
            cell_array[2, 2] = float(line.split()[2])
        elif i == 5:
            list_of_elements = line.split()
        elif i == 6:
            number_of_atoms = line.split()
            number_of_atoms = [int(x) for x in number_of_atoms]
            atom_positions_list = []
        elif 7 < i <= sum(number_of_atoms) + 7:
            line = [float(x) for x in line.split()]
            atom_positions_list.append(line)
        elif 7 < i > sum(number_of_atoms) + 7:
            break

# header to savetxt function
with open(poscar_file, "r") as pos_file:
    head = list(islice(pos_file, 8))

# array with atom positions in base structure
atom_positions = np.asarray(atom_positions_list)

# array with sqrt of masses for later calculations
masses_sqrt_array = np.full((3, number_of_atoms[0]), math.sqrt(element(list_of_elements[0]).mass))
for i in range(len(list_of_elements)-1):
    temporary_array = np.full((3, number_of_atoms[i+1]), math.sqrt(element(list_of_elements[i+1]).mass))
    masses_sqrt_array = np.concatenate((masses_sqrt_array, temporary_array), axis=1)

if create_movie is True:
    # lattice parameters
    Cell.a = np.linalg.norm(cell_array[0, :])
    Cell.b = np.linalg.norm(cell_array[1, :])
    Cell.c = np.linalg.norm(cell_array[2, :])
    Cell.alpha = angle_between(cell_array[1, :], cell_array[2, :])
    Cell.beta = angle_between(cell_array[0, :], cell_array[2, :])
    Cell.gamma = angle_between(cell_array[0, :], cell_array[1, :])

    # list of molecules for pdb file
    molecules_list = []
    list_index = 0
    for i in number_of_atoms:
        for j in range(number_of_atoms[list_index]):
            molecules_list.append(list_of_elements[list_index])
        list_index += 1

lines_with_freq = []

outcar_positions = np.zeros((6, sum(number_of_atoms)))

# reading outcar file
with open(outcar_file) as out_file:
    for line in out_file:
        if line.strip() == 'Eigenvectors and eigenvalues of the dynamical matrix':
            break
    for _ in range(3):
        next(out_file)
    for line in out_file:
        if line.strip() == '-----------------------------------------------------------------' \
                           '---------------------------------------':
            break
        if 'cm-1' in line:
            next(out_file)
            lines_with_freq.append(line)
            temporary_list = []
            temporary_array = np.zeros((6, sum(number_of_atoms)))
            for i in range(sum(number_of_atoms)):
                temporary_list.append(next(out_file).split())
            k = int(line.split()[0]) - 1
            for j in range(sum(number_of_atoms)):
                for i in range(6):
                    temporary_array[i, j] = temporary_list[j][i]
            if k != 0:
                outcar_positions = np.dstack((outcar_positions, temporary_array))
            else:
                outcar_positions = temporary_array

# limitation to shifts
outcar_positions = outcar_positions[3:6, :, :]

# list of desired modes
modes = []
for line in lines_with_freq:
    ind = (line.split()).index('cm-1')
    if not freqs:
        if freq_min <= float(line.split()[ind - 1]) <= freq_max:
            modes.append(line)
    else:
        if str("%.2f" % float(line.split()[ind - 1])) in freqs:
            modes.append(line)

# array of amplitude
amp_range = list('{:.2f}'.format(i) for i in np.arange(amp_min, amp_max + amp_step, amp_step))
amp_range = list(float(i) for i in amp_range)

deformed_array = np.copy(atom_positions)

# main calculations
for mode in modes:
    shifts_array = np.copy(outcar_positions[:, :, int(mode.split()[0])-1])
    divided_array = np.copy(shifts_array)
    np.divide(shifts_array, masses_sqrt_array, out=divided_array)
    ind = (mode.split()).index('cm-1')

    for amp in amp_range:
        normalised_with_amp = divided_array * amp
        normalised_with_amp = normalised_with_amp.T
        np.add(atom_positions, normalised_with_amp, out=deformed_array)

        # fast solution not to overwrite a file
        name_of_file = file_name + '_' + str(mode.split()[ind-1]) + '_0_' + str(amp) + '.vasp'
        if os.path.isfile(name_of_file) is True:
            name_of_file = file_name + '_' + str(mode.split()[ind-1]) + '_1_' + str(amp) + '.vasp'
            if os.path.isfile(name_of_file) is True:
                name_of_file = file_name + '_' + str(mode.split()[ind-1]) + '_2_' + str(amp) + '.vasp'
        header = ''.join(map(str, head))
        header = header.rstrip('\n')
        np.savetxt(name_of_file, deformed_array, delimiter=' ', fmt='%.15f', header=header, comments='')

        # movies
        if create_movie is True:
            position_parameters = -1.000, -0.9778, -0.9157, -0.8230, -0.7088, -0.5802, -0.4421, -0.2978, -0.1497, \
                                  0.0000, 0.1497, 0.2978, 0.4421, 0.5802, 0.7088, 0.8230, 0.9157, 0.9778, 1.0000
            movie_file_name = 'Movie_' + file_name + '_' + str(mode.split()[ind - 1]) + '_' + str(amp) + ".pdb"
            f = open(movie_file_name, 'w')
            model_numb = 1

            for par in position_parameters:
                parametrised_array = normalised_with_amp * par
                np.add(atom_positions, parametrised_array, out=deformed_array)
                np.around(deformed_array, decimals=3, out=deformed_array)
                f.write('MODEL' + "{:>9}".format(model_numb) + "\n")
                cryst_line = "CRYST1" \
                             + "{:>9.3f}".format(Cell.a) \
                             + "{:>9.3f}".format(Cell.b) \
                             + "{:>9.3f}".format(Cell.c) \
                             + "{:>7.2f}".format(Cell.alpha) \
                             + "{:>7.2f}".format(Cell.beta) \
                             + "{:>7.2f}".format(Cell.gamma)
                f.write(cryst_line + "\n")
                for i in range(sum(number_of_atoms)):
                    coord = np.zeros(3)
                    np.dot(cell_array, deformed_array[i, :], out=coord)
                    np.around(coord, decimals=3, out=coord)
                    atom_line = "ATOM" \
                                + "{:>7}".format(i+1) \
                                + "  " \
                                + "{:<2}".format(molecules_list[i]) \
                                + "   MOL          " \
                                + "{:>7}".format(coord[0]) \
                                + "{:>7}".format(coord[1]) \
                                + "{:>7}".format(coord[2]) \
                                + "  1.00  0.00          " \
                                + "{:<2}".format(molecules_list[i])
                    f.write(atom_line + "\n")
                f.write("ENDMDL\n")
                model_numb += 1
