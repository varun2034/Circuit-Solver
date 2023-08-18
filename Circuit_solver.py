"""
EE2703 Applied Programming Lab - 2022
Assignment 2: Solution
Jarpla Yashwanth - EE20B048
26th January 2021
"""

#INPUT: .netlist file
#OUTPUT: Identifies errors in SPICE program code, and displays the tokens in reverse order


from sys import argv, exit
from xml.dom.minicompat import NodeList
from numpy import *

if len(argv) == 2 :	
	file = argv[1]  #taking .netlist file in file variable
else:
	print("Number of arguments must be 2")  #gives the warning if file is not in required format
	exit(0)

# File Handling    
try:
	f = open(file)
	lines = f.readlines()
	f.close()
except:
	print("File not found")
	exit(0)

start_Cir = '.circuit'   # variable declaration for some strings in .netlist file
end_Cir = '.end'
AC = '.ac' 
start = -1              #initialising the variables start and END
end = -2
frequency = 0

# Now we define a class of all elements which can be used when they are called again in program.
class ALL_components():                            # General class for all the elements
    def __init__(self,line):
        self.line = line
        self.tokens = self.line.split()
        self.from_node = self.tokens[1]  #n1 and n2 values in the circuit are passed to the variables
        self.to_node = self.tokens[2]
        if self.tokens[0][0] ==  'H' or self.tokens[0][0] =='F':   #checking for controlled sources 
            self.type = 'currentControlled'
            self.value = float(self.token[4]) 
            self.vol_source = self.tokens[3]
        elif len(self.tokens) == 5:
            self.type = 'dc'                      #checking arguments for dc type circuit
            self.value = float(self.tokens[4])
        elif self.tokens[0][0] ==  'E' or self.tokens[0][0] =='G': 
            self.type = 'voltageControlled'
            self.value = float(self.tokens[5])
            self.vol_node1 = self.tokens[3]
            self.vol_node2 =self.tokens[4]
        elif len(self.tokens) == 6:
            self.type = 'ac'
            Vac = float(self.tokens[4])/2   #checking and analyising for AC type circuit
            phase = float(self.tokens[5])
            Real_Part = Vac*cos(phase)
            Imag_Part = Vac*sin(phase)
            self.value = complex(Real_Part,Imag_Part)
        else:
            self.type = 'RLC'                 #RLC circuit
            self.value = float(self.tokens[3])

def convert_to_dict(Token_line):
    ''' Returns a dictionary of nodes from the circuit definition. By default, the GND node is assigned value 0. '''
    dictionaryOfNodes = {}
    nodes = [ALL_components(line).from_node for line in Token_line]      		# adds all line's node1
    nodes.extend([ALL_components(line).to_node for line in Token_line])  		# adds all line's node2
    ind = 1
    nodes = list(set(nodes))					# using set remove all repeated nodes
    nodes = sorted(nodes)  
    for node in nodes:							# make  dictionary of all nodes as key and value is index
        if node == 'GND' :						# GND - index = 0
            dictionaryOfNodes[node] = 0
        else :
            dictionaryOfNodes[node] = ind
            ind += 1
    return dictionaryOfNodes


def insert_key(diction,value):
    ''' Gets the corresponding key for a value in the dictionary '''
    for key in diction.keys():
        if diction[key] == value :              #defining keys for nodes
            return key


def mkDictionary(Token_line,element):
    ''' Makes a dictionary for each component of the particular type of element ''' 
    e = element
    ele_dict = {}
    ele_names = [ALL_components(line).tokens[0] for line in Token_line if ALL_components(line).tokens[0][0].lower()== e]
    ind=0
    for name in ele_names:
        ele_dict[name] = ind      #defining value for a key
        ind += 1
    return ele_dict

def FinalMatrix(M,b,Token_line,volt_dict,vcvs_dic,ccvs_dic,dictionaryOfNodes):
    for ele in Token_line:
        idx = Token_line.index(ele)
        element = ALL_components(Token_line[idx])
        element_name = Token_line[idx].split()[0]
        nodeL = dictionaryOfNodes[element.from_node]    #1st node
        nodeK = dictionaryOfNodes[element.to_node]      #2nd node
        value = element.value
        if element_name[0] == 'R':
            M[nodeK][nodeK]+= 1/value
            M[nodeL][nodeL]+= 1/value
            M[nodeL][nodeK]-= 1/value
            M[nodeK][nodeL]-= 1/value
        if element_name[0] == 'L':
            M[nodeK][nodeK]-= complex(0,1/(2*pi*frequency*value))
            M[nodeL][nodeL]-= complex(0,1/(2*pi*frequency*value))
            M[nodeL][nodeK]+= complex(0,1/(2*pi*frequency*value))
            M[nodeK][nodeL]+= complex(0,1/(2*pi*frequency*value))
        if element_name[0] == 'C':
            M[nodeK][nodeK]+= complex(0, 2*pi*frequency*value)
            M[nodeL][nodeL]+= complex(0, 2*pi*frequency*value)
            M[nodeL][nodeK]-= complex(0, 2*pi*frequency*value)
            M[nodeK][nodeL]-= complex(0, 2*pi*frequency*value)
        if element_name[0] == 'G':  #vccs
            nodeM = element.vol_node1
            nodeN = element.vol_node2
            M[nodeK][nodeM]+=value
            M[nodeK][nodeN]-=value
            M[nodeL][nodeM]-=value
            M[nodeL][nodeN]+=value
        if element_name[0] == 'H': #ccvs
            nodeKL = len(dictionaryOfNodes)+len(volt_dict)+len(vcvs_dic)+ccvs_dic[element.tokens[0]]
            M[nodeK][nodeKL]+=1
            M[nodeK][nodeKL]-=1
            M[nodeKL][nodeK]-=1
            M[nodeKL][nodeL]-=1
            M[nodeKL][len(dictionaryOfNodes)+volt_dict[element.vol_source]]+=value
        if element_name[0] == 'E':  #vcvs
            nodeM = element.vol_node1
            nodeN = element.vol_node2
            nodeKL = len(dictionaryOfNodes)+len(volt_dict)+vcvs_dic[element.tokens[0]]
            M[nodeK][nodeKL]+=1
            M[nodeL][nodeKL]-=1
            M[nodeKL][nodeK]+=1
            M[nodeKL][nodeL]-=1
            M[nodeKL][nodeM]-=value
            M[nodeKL][nodeN]=value

        if element_name[0] == 'F':  #cccs
            nodeMN = len(dictionaryOfNodes)+ volt_dict[element.vol_source]
            M[nodeK][nodeMN] += value
            M[nodeL][nodeMN]-= value
        if element_name[0] == 'V':
            nodeKL = len(dictionaryOfNodes) + volt_dict[element.tokens[0]]
            nodeK ,nodeL =nodeL,nodeK    
            M[nodeK][nodeKL]+=1
            M[nodeL][nodeKL]-=1
            M[nodeKL][nodeL]-=1
            M[nodeKL][nodeK]+=1
            b[nodeKL]+= value
        if element_name[0] == 'I':
            b[nodeK]+= value
            b[nodeL]-=value

    return M,b
	
    

for line in lines:

	location = lines.index(line)
	line = line.replace('\n','')
	line = line.split('#')[0]
	line = line.replace('\t',' ')	
	line=line.strip()
	lines[location] = line

	if line[:len(start_Cir)] == start_Cir:
		start = location
	elif line[:len(end_Cir)] == end_Cir:
		end = location;	
	elif line[:len(AC)] == '.ac' :
		frequency = float(line.split()[2])

if start >= end:
	print("Given Circuit Data Invalid")
	exit(0)


Token_line = [] 
for line in lines[start+1:end]:
	
	location = lines.index(line)
	line_list = line.split(" ")
	line_list = [elem for elem in line_list if elem != ""]

	if line_list == []:
		continue	

	if line_list[0][0] == 'R' or line_list[0][0] == 'L' or line_list[0][0] =='C' or line_list[0][0] =='I' :	
		if len(line_list) != 4:
			print("Incorrect Number of Parameters: Line ",location)
			exit(0)
		if line_list[1].isalnum() != True or line_list[2].isalnum() != True :
			print("Incorrect Node Representation - only alphanumeric variables: Line ",location)
			exit(0)
		
	elif line_list[0][0] ==  'E' or line_list[0][0] =='G':
		if len(line_list) != 6:
			print("Incorrect Number of Parameters: Line ",location)
			exit(0)
		if line_list[1].isalnum() != True or line_list[2].isalnum() != True or line_list[3].isalnum() != True or line_list[4].isalnum() != True:
			print("Incorrect Node Representation - only alphanumeric variables: Line ",location)
			exit(0)
	

	elif line_list[0][0] ==  'H' or line_list[0][0] =='F':
		if len(line_list) != 5:
			print("Incorrect Number of Parameters: Line ",location)
			exit(0)
		if line_list[1].isalnum() != True or line_list[2].isalnum() != True:
			print("Incorrect Node Representation - only alphanumeric variables: Line ",location)
			exit(0)
		if line_list[3][0] != 'V':
			print("Incorrect Voltage Label: Line ",location)
			exit(0)

	Token_line.append(line)

dictionaryOfNodes = convert_to_dict(Token_line)
volt_dict = mkDictionary(Token_line,'v')
vcvs_dic = mkDictionary(Token_line,'e')
ccvs_dic = mkDictionary(Token_line,'h')
ind_dict = mkDictionary(Token_line,'l')

n = len(dictionaryOfNodes)
k = len(volt_dict)

dim = n + k + len(vcvs_dic) + len(ccvs_dic)  #dimension of matrix : length of node + length of voltage source(indep. +vcvs + ccvs)

M = zeros((dim,dim), dtype=complex)
b = zeros(dim, dtype=complex)
 
if frequency == 0:     # if circuit is dc
	M = zeros((dim+len(ind_dict),dim+len(ind_dict)),dtype=float)
	b = zeros(dim+len(ind_dict),dtype=float)

M,b = FinalMatrix(M,b,Token_line,volt_dict,vcvs_dic,ccvs_dic,dictionaryOfNodes)

M[0] = 0                # make first row zero
M[0,0] =1               
	
try:
	x = linalg.solve(M,b)    
except Exception:
	print('The incidence matrix cannot be inverted as it is singular. Please provide a valid circuit definition')
	exit()
# display values
for i in range(len(dictionaryOfNodes)):
    print("The voltage at node {} is {}".format(insert_key(dictionaryOfNodes,i),x[i]))
for j in range(len(volt_dict)):
    print('The current through source {} is {}'.format(insert_key(volt_dict,j),x[n+j]))
for p in range(len(vcvs_dic)):
    print('The current through source {} is {}'.format(insert_key(vcvs_dic,p),x[n+k+p]))
for q in range(len(ccvs_dic)):
    print('The current through source {} is {}'.format(insert_key(ccvs_dic,q),x[n+k+q+len(vcvs_dic)])) 
	

