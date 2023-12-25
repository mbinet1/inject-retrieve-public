### All constants are defined here

# astrophysical constants
G = 6.673e-11  # gravitational constant [m^3 / kg / s^2]
R_sun = 695508000  # radius of the Sun [m]
R_earth = 6371000  # radius of the Earth [m]
earth_sun_ratio = 0.00916
R_jup = 69911000  # radius of Jupiter [m]
M_sun = 1.989 * 10 ** 30  # mass of the Sun [kg]
seconds_per_day = 86400

earth_transit_dur = 0.52
days_in_year = 365.25

kepler_Q2_start = 170
kepler_Q3_end = 345.573
kepler_Q2_end = 255.573

plato_time_origin = 171072143.75
plato_points_per_quarter = 304128
plato_Q23_end, plato_Q24_start, plato_Q24_end = 88, 90, 178

raw_plato_cadence = 25

inf = float('inf')

# IDs of PLATO stars that are less than 5% off both the Sun's mass and radius
plato_sunlike_stars = [9, 71, 106, 173, 190, 271, 282, 287, 329, 330, 349, 381, 396, 465, 482, 537, 539, 543, 557, 578, 601, 624, 646, 664, 669, 680, 713, 740, 790, 832, 900, 955, 974]
# V > 10
plato_sunlike_faint = [71, 329, 330, 349, 465, 482, 537, 557, 578, 624, 646, 664, 669, 680, 740, 790, 832]
# V < 10
plato_sunlike_bright = [9, 106, 173, 190, 271, 282, 287, 381, 396, 539, 543, 601, 713, 900, 955, 974]

# V > 10
kepler_sunlike = {11: [4903847, 11138101, 11128041, 6063284, 3733537, 9356691, 8408982, 8145383, 5788360],
                        10: [8038388, 5357446, 4243697, 8424992, 9898385, 9468847],
                        9: [7801848, 3427720]}


#kepler nezw mag10 = [4554810, 4255487, 9693187, 11029516, 4921578, 7090654, ]
#kepler_not_sunlike = [9468847, 8145383, 5788360, 5949536, 9098294]
#kepler_not_sunlike_at_all = [7335019, 9307354, 10934586, 9764820, 9788260, 8216936]

all_plato_cameras = range(1,25)

plato_total_nb_stars = 1000

# A few stars of the plato_sunlike_stars selection, with informations for one or two quarters
# Format: {star_id_quarter: [star_id, quarter, rounded apparent mag, number of cameras, noise in ppm]}:
star_params = {"832_Q24": [832, "Q24", 11, 8, 859], 
        "646_Q23": [646, "Q23", 11, 14, 738],
        "646_Q24": [646, "Q24", 11, 24, 581],
        "329_Q23": [329, "Q23", 10, 6, 883],
        "329_Q24": [329, "Q24", 10, 12, 647],
        "624_Q23": [624, "Q23", 10, 24, 418],
        "624_both" : [624, "both", 10, 24, 418],
        "790_both" : [790, "both", 10, 24, 450],
        "900_Q24": [900, "Q24", 9, 6, 523], 
        "9_Q23": [9, "Q23", 9, 12, 354],
        "396_Q23": [396, "Q23", 9, 24, 241]}

## Noise of stars!!
## Number of cams per star

#Star number 9: [12, 6]
#Star number 71: [7, 12]
#Star number 106: [12, 12] 9.05
#Star number 173: [4, 12]
#Star number 190: [12, 12] 8.98
#Star number 271: [11, 18]
#Star number 282: [12, 6]
#Star number 287: [12, 12] 9.11
#Star number 329: [6, 12]
#Star number 330: [22, 24] 10.09
#Star number 349: [6, 24]
#Star number 381: [12, 24]
#Star number 396: [24, 24] 8.95
#Star number 465: [11, 24]
#Star number 482: [24, 18]
#Star number 537: [18, 16] 10.68
#Star number 539: [6, 24]
#Star number 543: [12, 24]
#Star number 557: [12, 18]
#Star number 578: [24, 21] 10.37
#Star number 601: [12, 18]
#Star number 624: [24, 24] 10.02
#Star number 646: [14, 23]
#Star number 664: [24, 6]
#Star number 669: [12, 19]
#Star number 680: [16, 12]
#Star number 713: [24, 24] 9.41
#Star number 740: [17, 12]
#Star number 790: [22, 24] 10.39
#Star number 832: [10, 8] 
#Star number 900: [0, 6]
#Star number 955: [0, 6]
#Star number 974: [0, 6]


