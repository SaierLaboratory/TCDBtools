'''
 * FileName: formatQuodDomain.py
 * Author: Yichi Zhang
 * Description: this file is a module for quod.py in Saier's lab. The module will contain a main function
              and several private functions to format the domain part string in quod.py(start from -ar), it suppose to
              set layers of protein domains and colors, and return a substring for drawing domains in quod.py
 * Input: qacc/sacc - query/subject accessions
	  dic: *************************************IMPORTANT******************************************
		main functions of the module only accept dictionaries in a certain format
		dictionary D is a dictionary containing accessions as keys, and lists of lists as values. 
		Example:
			


                                                            
                              ____protein accession 1 (key): list of lists L1 
                              | 
                              |___protein accession 2      : list of lists L2
			D_____|
			      |___protein accession 3      : list of lists L3
				...


			L-----[ sublist l1, sublist l2, sublist l3 ... ]
			
			Each sublist l contains all necessary information about a domain of that protein
			
			l-----[ envlope_from(int), envlope_to(int), domain_acc(str), clan_accession(str), 0(color index), 0(layer index)........]
				|__________________________________________________________________________________________|
								             |
				The first 6 elements of l must be in this format. Users are welcome to add more information
				to the dictionary, elements after the 4th are not relavent to module functions. If the domain
				has no clan_acc, then the 4th element must be 'N/A'

		************************************************************************************************
	     bottom: specify the postion of the lowest domain on the plot
	     margin: specify the gap between 2 domains on the plot
 * Output: 1. a string for quod.py -ar option (single protein)
	   2. a list constains 2 strings for quod.py -ar option( protein comparision)
'''
class Domain_string(object):
	__colors=['red','teal','blue','yellow','purple','maroon','black','pink','green','cyan','beige','lavendar','orange','grey',
	  'olive','navy','magenta','coral','lime','plum','chocolate','indigo','brown']
	def __init__(self):
		__bottom = 0
		__margin = 0

	def __setBottom( self, number ):
		self.__bottom = number

	def __setMargin( self, number ):
		self. __margin = number

	def __setLayer( self, plist ):
        # in a nested loop, check if domains are overlap. overlapped domains will be elevated to a higher layer
		i = 0
        	while ( i < len(plist)-1 ):
            		k = i + 1
            		while ( k < len(plist) ):
                		if int(plist[k][0]) < int(plist[i][1]):
                		    	if plist[k][5] == plist[i][5]:
                	        		plist[k][5] = plist[k][5] + 1
                    		k = k + 1

            		i = i + 1

    	def __setColor_double( self, plist1, plist2 ): #plist1 is the pfam list for query, plist2 is the pfam_list for subject
        # pick all clans from the query and subject. Once a unique clan is encountered, store it into the list. the index of that clan will also be the index of the color in the color list, and it will be appended to the end of the pfam_dic lists
        	clan_list = []
        	for x in plist1:
			if x[3] != 'N/A':
        	        	if x[3] not in clan_list:
        	                	clan_list.append(x[3])
				x[4] = clan_list.index(x[3])
			else:
				if x[2] not in clan_list:
					clan_list.append(x[2])
				x[4] = clan_list.index(x[2])
		for y in plist2:
			if y[3] != 'N/A':
        	        	if y[3] not in clan_list:
        	                	clan_list.append(y[3])
				y[4] = clan_list.index(y[3])
			else:
				if y[2] not in clan_list:
					clan_list.append(y[2])
				y[4] = clan_list.index(y[2])

	def __setColor_single( self, plist ): # fucntion overloading. It draws colors for a single protein
	        clan_list = []
        	for x in plist:
			if x[3] != 'N/A':
        	        	if x[3] not in clan_list:
        	                	clan_list.append(x[3])
				x[4] = clan_list.index(x[3])
			else:
				if x[2] not in clan_list:
					clan_list.append(x[2])
				x[4] = clan_list.index(x[2])

	def two_proteins( self, qacc, sacc, pfam_dic, bottom = -2.8, margin = 0.3 ):# for comparing two proteins together
		self.__setBottom( bottom )
        	self.__setMargin( margin )
        	result = []
	# sort qacc list and sacc list in the dictionary in accending form based on the left coordinates
	# coordinates are string numbers
		pfam_dic[qacc].sort( key=lambda x: int(x[0]) )
		pfam_dic[sacc].sort( key=lambda x: int(x[0]) )
	# in each list, decide whether domains are overlapped and whether they have the same color to paint append their colors and layers to the end of the sublist
		self.__setColor_double( pfam_dic[qacc], pfam_dic[sacc] )
		self.__setLayer( pfam_dic[qacc] )
		self.__setLayer( pfam_dic[sacc] )
	# try to produce domain strings and concanate to the end of the result
		q_str = '-ar '
		s_str = '-ar '
		for x in pfam_dic[qacc]:
			color = self.__colors[x[4]]
			q_str = q_str + '{f}-{t}:"{text}":{b}:{c} '.format( f=x[0],t=x[1],text=x[2],b=self.__bottom+self.__margin*x[5], c=color )
		for y in pfam_dic[sacc]:
			color = self.__colors[y[4]]
			s_str = s_str + '{f}-{t}:"{text}":{b}:{c} '.format( f=y[0],t=y[1],text=y[2],b=self.__bottom+self.__margin*y[5], c=color )

		result.append(q_str)
		result.append(s_str)
		return result

	def one_protein( self, qacc, pfam_dic, bottom = -2.8, margin = 0.3 ):# function overloading, for single protein plot
        	self.__setBottom( bottom )
        	self.__setMargin( margin )
        	# sort qacc list and sacc list in the dictionary in accending form based on the left coordinates
        	# coordinates are string numbers
        	pfam_dic[qacc].sort( key=lambda x: int(x[0]) )
        	# assign colors according to number of different clans
        	self.__setColor_single( pfam_dic[qacc] )
		self.__setLayer( pfam_dic[qacc] )
        	# try to produce domain strings and concanate to the end of the result
        	q_str = '-ar '
        	for x in pfam_dic[qacc]:
        	        color = self.__colors[x[4]]
        	        q_str = q_str + '{f}-{t}:"{text}":{b}:{c} '.format( f=x[0],t=x[1],text=x[2],b=self.__bottom+self.__margin*x[5], c=color )
        	return q_str

