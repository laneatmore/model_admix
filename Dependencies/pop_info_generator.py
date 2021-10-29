#!/usr/bin/env python3


#generates pop id's

import subprocess
import sys
import os
import fileinput

os.chdir("../Admixture/")

N1 = int(sys.argv[1])
N2 = int(sys.argv[2])
N3 = int(sys.argv[3])

def new_file(N1, N2, N3):
	if (len(sys.argv) < 4):
		print("not enough arguments")
		exit()
	elif (len(sys.argv) > 5):
		print("too many arguments")
		exit()
	else:
		pass
	
	with open("population_information.txt", "w") as file: 
		file.write(str("ID   col2    col3     POP     REGION\n"))
		for i in range(N1):
			my_id = "msp_" + str(i)
			file.write(my_id + str("  0     PAR1     0     PAR1\n"))
		for i in range(N2):
			my_id2 = "msp_" + str(N1 + i)
			file.write(my_id2 + str("  1     PAR2     1     PAR2\n"))	
		for i in range(N3):
			my_id3 = "msp_" + str(N1 + N2 + i)
			file.write(my_id3 + str("  2      ADM       2     ADM\n"))
	
	with fileinput.FileInput("population_information.txt", inplace = True,
		backup=".bak") as file:
			for line in file:
				print(line.replace("msp_0", "msp_00"), end='')
		
					
new_file(N1, N2, N3)
		