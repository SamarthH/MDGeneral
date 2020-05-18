ip_file = open("input.txt", "r")
op_file = open("input_params.txt", "w")

con = ip_file.readlines()
con = [i.strip('\n') for i in con]

for index, i in con:
	if(index%2 == 0):
		continue
