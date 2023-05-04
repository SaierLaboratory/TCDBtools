import numpy as np

def collapse(arr, nanlike=np.nan):
	''' eliminate all nans, spitting out the original indices '''
	if isinstance(nanlike, float) and np.isnan(nanlike): 
		newarr = arr[~np.isnan(arr)]
		nans = np.nonzero(np.isnan(arr))[0]
	else: 
		arr = np.array(arr)
		newarr = arr[arr != nanlike]
		nans = np.nonzero(arr == nanlike)[0]

	return newarr, nans

def uncollapse(arr, nans, nanlike=np.nan, ignore_oob=False):
	''' add back all nans given indices '''
	if not len(nans): return arr
	if max(nans) > (len(arr) + len(nans)) and not ignore_oob: raise IndexError('Found a nan with an index outside the allowed range')
	nanset = set(nans)

	newarr = np.zeros(len(arr) + len(nans))
	j = 0
	for i in range(len(newarr)):
		if i in nanset: 
			try: newarr[i] = nanlike
			except IndexError as e:
				if ignore_oob: continue
				else: raise e
		else:
			try: newarr[i] = arr[j]
			except IndexError as e:
				if ignore_oob: continue
				else: raise e
			j += 1

	return newarr

def rle(arr):
	pass
