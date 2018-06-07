from sys import argv
import os

script, filename = argv

def replace_number_of_bonds():
    with open(filename, 'r') as input_file, open('a.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '26 26  0  0  0  0            999 V2000':
                output_file.write(' 26 24  0  0  0  0            999 V2000\n')
            else:
                output_file.write(line)

def replace_1st_bond():
    with open('a.sdf', 'r') as input_file, open('b.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '2  3  2  0  0  0  0':
                output_file.write('  2  3  1  0  0  0  0\n')
            else:
                output_file.write(line)

def replace_2nd_bond():
    with open('b.sdf', 'r') as input_file, open('c.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '3  4  1  0  0  0  0':
                output_file.write('  3  4  2  0  0  0  0\n')
            else:
                output_file.write(line)

def replace_3rd_bond():
    with open('c.sdf', 'r') as input_file, open('d.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '5  6  1  0  0  0  0':
                output_file.write('  5  6  2  0  0  0  0\n')
            else:
                output_file.write(line)

def replace_4th_bond():
    with open('d.sdf', 'r') as input_file, open('e.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '1  2  1  0  0  0  0':
                output_file.write('  1  2  2  0  0  0  0\n')
            else:
                output_file.write(line)

def delete_5th_bond():
    with open('e.sdf', 'r') as input_file, open('f.sdf', 'w') as output_file:
        for line in input_file:
            if line.strip() == '4  5  1  0  0  0  0':
                output_file.write('')
            else:
                output_file.write(line)

def delete_6th_bond():
    with open('f.sdf', 'r') as input_file, open(filename, 'w') as output_file:
        for line in input_file:
            if line.strip() == '1  6  1  0  0  0  0':
                output_file.write('')
            else:
                output_file.write(line)

replace_number_of_bonds()

replace_1st_bond()

replace_2nd_bond()

replace_3rd_bond()

replace_4th_bond()

delete_5th_bond()

delete_6th_bond()

os.remove('a.sdf')
os.remove('b.sdf')
os.remove('c.sdf')
os.remove('d.sdf')
os.remove('e.sdf')
os.remove('f.sdf')





