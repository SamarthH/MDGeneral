ip_file = open("input.txt", 'r')
input = ip_file.read()

output = ""

arr = input.split('\n')
letters = [10, 15, 14, 11, 10, 14, 7, 23, 20]

for index, value in enumerate(arr):
  output+=value[letters[index]:]
  output+=","

op_file = open("pro_input.txt", "w")
op_file.write(output)

op_file.close()
ip_file.close()
