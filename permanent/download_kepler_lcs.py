import dill 
import pickle
import numpy as np 
import csv
import lightkurve as lk

quarters = [2,3]

class Data:
	pass

data = Data()
delta = 2454833
k = 0
star_ids = [9307354, 4903847, 10934586, 11395018, 4921578, 
            7500161, 11138101, 11182897, 11128041, 5357446, 
            4753390, 10160534, 3733537, 4243697, 9764820, 
            8408982, 8424992, 9788260, 8216936, 9898385, 
            4660971, 9468847, 8145383, 5788360, 8038388,
            7335019, 10908196]
data.stars_used = star_ids
nb_stars_total = len(star_ids)
indices_kept = range(3)
nb_stars_kept = len(indices_kept)

def GetMask(curve, i):
        mask1 = (curve.time.value>0)
        mask2 = mask1
        if i==0: 
            start, end = 700, 850
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==1:
            start, end = 100, 180
            start1, end1 = 575, 595
            start2, end2 = 1100, 1130
            return (curve.time.value<start) | (curve.time.value>end), (curve.time.value>end1) | (curve.time.value<start1), ((curve.time.value>end2) | (curve.time.value<start2)), True
        elif i==2:
            start, end = 700, 770
            start1, end1 = 1050, 1130
            return (curve.time.value < start) | (curve.time.value > end), (curve.time.value>end1) | (curve.time.value<start1), mask2, True
        elif i==4:
            start, end = 290, 360
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==5:
            start, end = 750, 950
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==8:
            start, end = 1150, 1200
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2,True
        elif i==9: 
            start, end = 720, 980 ###New try
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2,True
        elif i==11:
            start, end = 300, 350
            start1, end1 = 1100, 1140
            return (curve.time.value<start) | (curve.time.value>end), (curve.time.value>end1) | (curve.time.value<start1), mask2, True
        elif i==12: 
            start, end = 0, 200 ###New try
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==14: 
            start, end = 1200, 1280
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==15: 
            start, end = 1180, 1200
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==17:
            start, end = 700, 750
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        elif i==19:
            start, end = 230, 270
            return (curve.time.value < start) | (curve.time.value > end), mask1, mask2, True
        return None, None, None, False

data.all_lc = []
data.all_timestamps = []

for i in indices_kept:
    data.all_timestamps.append([])
    star_id = star_ids[i]
    print("Loading star n°", i, " : KIC",star_id)
    search_result = lk.search_lightcurve("KIC"+str(star_id), author = 'Kepler', cadence = 'long', radius=0, quarter=quarters)
    lc_collection = search_result.download_all()
    new_collection = []
    for j in range(len(lc_collection)):
        mask, mask1, mask2, exists = GetMask(lc_collection[j], i)
        if exists:
            lc_collection[j] = lc_collection[j][mask]
        nb_points = len(lc_collection[j])
        if nb_points>0:
            first = lc_collection[j].time[0].to_value('jd','float') - delta
            last = lc_collection[j].time[-1].to_value('jd','float') - delta
            data.all_timestamps[k].append([first,last,nb_points])
    lc = lc_collection.stitch()
    #lc.plot();
    data.all_lc.append(lc)
    k+=1
print("Loading light curves finished")



with open('all_kepler_lightcurves_Q2-Q3.pkl', 'wb') as f:
    pickle.dump(data, f)
    f.close
