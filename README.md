# phonons.py


A script that reads POSCAR and OUTCAR files from VASP package and produces VASP files with deformed structures for specific frequency and amplitude of vibration. 

## installation
Script can be launched from any directory if it is installed in such a way:
1. Download phonons.py file
2. Put it in a folder in your $PATH (such as ```usr/local/bin```)
3. Change permissions by typing ```chmod +x phonons.py; mv phonons.py phonons```

Other option is to have the script in your working directory and launching it using command ```python phonons.py```.

## basic usage

Script can be operated in manual mode or using input file with specified parameters. The second option is recommended for generating structures for many frequencies in one run. This mode starts automatically if the run command is followed by an input file name. Example:
```
phonons input.txt
```

#### input file compositon
An example file is available in this repository. Comment lines with \# signs must remain in place. The file layout is as follows:
```
1. comment line
2. framework name - it is name of POSCAR and OUTCAR file you are using (without extension), both files must be named the same
3. comment line
4. desired amplitude range (start stop step), e.g., -2 2 0.5 (spaces as separators),  start and stop values are included in interval
5. comment line
6. creation of a PDB file with vibration movement, True or False (False is default which means that PDB file is not created)
7. comment line
8. list of desired frequencies (in cm-1) separated with space and rounded up to two decimals 
```

#### manual mode
If no input file is specified user will be asked to enter the values manually. 
To exit program user must press enter without entering any values. 
