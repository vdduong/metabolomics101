####### HERE BEGIN the routines for k-means algorithm :

def list_point_append(dict_assignment, peak_list, pattern_tocsy) :
  """
	create the list of secondary chemical shifts
	scattering (X,Y) in order to take the k-mean algo in count"""
	list_point = []
	for key_assignment in dict_assignment.keys():
		peak_tocsy = pattern_tocsy[key_assignment] # name the peak in tocsy
		cs_tocsy_1 = float(peak_tocsy[0]) # assigned chemical shifts
		cs_tocsy_2 = float(peak_tocsy[1])
		item_assignment = dict_assignment[key_assignment]
		if len(item_assignment) > 0 :
			for item in item_assignment :
				cs_exp_1 = float(peak_list[item[0]][0])
				cs_exp_2 = float(peak_list[item[0]][1])
				coord_1 = cs_tocsy_1 - cs_exp_1
				coord_2 = cs_tocsy_2 - cs_exp_2
				coord_1 = float(('%.4f')%(coord_1))
				coord_2 = float(('%.4f')%(coord_2))
				list_point.append((coord_1,coord_2,key_assignment,item[0]))
				### scs_1, scs_2, position in pattern tocsy, position in peak_list
		else : 
			pass
	return list_point

def distanceBetweenPoints(peak_1, peak_2, tol_H):
	distance = math.sqrt((peak_1[0] - peak_2[0])**2 + \
					(peak_1[1] - peak_2[1])**2)/ tol_H
	return distance
	
def assignment_centre(list_centre, list_point, tol_H) :
	"""
	classify the elements of list_point into a dictionary
	whose keys are centres"""
	dict_centre = dict()
	for item in list_centre :
		dict_centre[item] = list() # initialization of dictionary of centres
	
	for peak_1 in list_point :
		min_local = 100.0
		listDistanceLocal = []
		for key_point in dict_centre.keys() :
			listDistanceLocal.append((distanceBetweenPoints(peak_1, key_point,tol_H), \
														key_point))
			if distanceBetweenPoints(peak_1, key_point, tol_H) < min_local :
				min_local = distanceBetweenPoints(peak_1, key_point, tol_H)
		for item in listDistanceLocal :
			if item[0] == min_local :
				dict_centre[item[1]].append(peak_1)
	return dict_centre
	
def update_centre(dict_centre) :
	"""update the centres until a stationary state"""
	list_centre = list()
	for key_centre in dict_centre.keys() :
		new_coord_1 = 0.0
		new_coord_2 = 0.0
		item_centre = dict_centre[key_centre] #it's a cluster of points
		for point in item_centre :
			new_coord_1+=point[0]
			new_coord_2+=point[1]
		if len(item_centre) > 0 :
			new_coord_1 = '%.4f'%(new_coord_1/float(len(item_centre)))
			new_coord_2 = '%.4f'%(new_coord_2/float(len(item_centre)))
		new_coord_1 = float(new_coord_1)
		new_coord_2 = float(new_coord_2)
		new_key_centre = (new_coord_1, new_coord_2)
		list_centre.append(new_key_centre)
	return list_centre	
	
def clustering_assignment(list_point, tol_H) :
	"""
	randomly chose the centres from list_point
	to fit the number k in k_means
	then 
	make the k_means algorithm
	assignment_centre
	update_centre
	"""
	nb_stab_limit = 100 # set a limit for ones that never converge
	filter_clustering = list() # to regroup the best clustering group of points
	filter_clustering_complet = list() # to regroup over the different clusters
                          # to form the best combinations
	k = round(math.sqrt(float(len(list_point))/2.0)) # thumb number formula
	
	list_centre = list()
	i = 1
	while i <=k :
          chococo = choice(list_point)
          if chococo not in list_centre :
            list_centre.append(chococo)
            i+=1
	nb_time_clustering = 1
	while 1 :
		dict_centre = assignment_centre(list_centre, list_point, tol_H)
		list_centre_1 = update_centre(dict_centre)
		if list_centre_1 == list_centre : 
			break
			#print 'situation stablized after %i times around %i centre'%(nb_time_clustering, \
											#len(list_centre))
		else : 
			list_centre = list_centre_1
			
		nb_time_clustering+=1
		if nb_time_clustering <= nb_stab_limit : pass
		else : 
			break
			#print 'situation has no convergent pattern'
	longest_element = 0
	for key in dict_centre.keys():
		if len(dict_centre[key]) > longest_element :
			longest_element = len(dict_centre[key])
	for key in dict_centre.keys():
		if len(dict_centre[key]) == longest_element :
			item_centre = dict_centre[key]
			for item_local in item_centre :
				filter_clustering.append((item_local[2], item_local[3]))
        for item in filter_clustering :
          filter_clustering_complet.append(item) # initialization by including all best elements
                                  # in the total group
        for key in dict_centre.keys():
          item_centre = dict_centre[key]
          for item_local in item_centre :
            if item_local[2] not in [item[0] for  item in filter_clustering] :
              filter_clustering_complet.append((item_local[2], item_local[3]))
        #print filter_clustering
        #print filter_clustering_complet
        #print '*'*5
        return filter_clustering, filter_clustering_complet
        
def distance_(list_1, list_2) :
	distance = 0.0
	for item_1 in list_1 :
		atom_local_1 = item_1[0]
		cs_local_1 = item_1[1]
		for item_2 in list_2 :
			if item_2[0] == atom_local_1 :
				cs_local_2 = item_2[1]
				if 'H' in atom_local_1 :
					distance+= (cs_local_1 - cs_local_2)**2/0.5**2
				elif 'C' in atom_local_1 :
					distance+= (cs_local_1 - cs_local_2)**2/0.98**2
				elif 'N' in atom_local_1 :
					distance+= (cs_local_1 - cs_local_2)**2/2.5**2
	distance = math.sqrt(distance)
	return distance
