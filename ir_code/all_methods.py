import current_parameters as prm
import numpy as np
import csv
import lightkurve as lk
import time
import pandas as pd
import pickle
import batman
import constants
import os
import shutil
if prm.is_notebook == False:
    if prm.transit_search == "gpu":
        from gpu_transit_toolkit import TCEDetector, TransitModeller
    elif prm.transit_search == "tls":
        from transitleastsquares import transitleastsquares


############################# CLASS AND FUNCTIONS !!!! ####################################################################################


class InjectRetrieve():
    
####################################### PART 1: functions that are only run once each time we start main program #######
    
    def __init__(self):
        np.random.seed(prm.numpy_seed)
        
        ## Computed amounts
        if prm.override_nb_transits:
            self.min_nb_transits, self.max_nb_transits = prm.min_nb_transits, prm.max_nb_transits
        else:
            if prm.transit_search == "gpu":
                self.min_nb_transits, self.max_nb_transits = 1, 1
            else: 
                self.min_nb_transits, self.max_nb_transits = 2, constants.inf

        if prm.mission == "plato":
            self.bin_factor = prm.plato_cadence//constants.raw_plato_cadence
            self.nb_stars_used = len(prm.plato_stars)
            self.initial_nb_points = constants.plato_points_per_quarter // self.bin_factor   
            if prm.bls_adapt_speed:
                prm.bls_frequency_factor *= 1/self.bin_factor   
        else:
            self.star_mag = prm.kepler_mag
            self.star_nb_cams = 1
            self.bin_factor = 1
            self.quarter_message = None
            self.chosen_quarter = 'both'
            self.nb_stars_used = len(prm.kepler_stars)
            self.initial_nb_points = constants.plato_points_per_quarter                                       
        
        ## Notebooks
        self.ntbk_str = ""
        if prm.is_notebook:
            self.ntbk_str = prm.notebook_string
        self.count_transits = 0
        # To compute then show the mean of correctly measured periods
        self.correct_average = 0   
        self.false_pos_average = 0   
        self.below_thresholds_avg = 0              
        self.sum_computing_transit_search = 0
        self.sum_computing_detrending = 0
        self.depth_ratios = []
        self.avg_depth_ratios = []
        self.dur_ratios = []

    def get_next_filename(self):
        counter = 1
        while True:
            filename = f""
            filename = f'../results/{self.base_filename}_{counter}.csv'
            if not os.path.exists(filename):
                return f"{self.base_filename}_{counter}"
            counter += 1
        return None

    def OpenWriter(self):
        if prm.base_filename == "default":
            if prm.mission == 'plato':
                self.base_filename = f'{prm.mission}_bin{prm.plato_cadence}'
                if prm.window_factor != 15:
                    self.base_filename += f'_wlf{prm.window_factor}'
                self.base_filename += f'_{prm.transit_search}'
                if self.nb_stars_used == 1:
                    self.base_filename += f'_star{prm.plato_stars[0]}'
            if prm.mission == 'kepler':
                self.base_filename = f'{prm.mission}'
                if prm.window_factor != 15:
                    self.base_filename += f'_wlf{prm.window_factor}'
            self.base_filename += f'_seed{prm.numpy_seed}'
            self.final_filename = self.get_next_filename()
        else:
            self.base_filename = prm.base_filename
        self.fgh = open(f'../results/{self.final_filename}.csv', 'w')
        self.writer = csv.writer(self.fgh)
        self.writer.writerow(["numpy_seed", "sample", "star_id", "quarter", "star_mag", "nb_cameras", 
                        "transit_search", "injected_radius", "impact_parameter", "transit_type"
                        "nb_transits_injected", 'nb_transits_recovered', "period_injected", "period_recovered", "relative_period_error", 
                        "t0_injected", "t0_recovered", "t0_error", "current_TPR"
                        'star_std', 'final_nb_points', 'new_SNR', "SNR", "search_power", "outcome",
                        "tdur_injected", "tdur_recovered", "tdur_relative_error",
                        "max_depth_injected", "mean_depth_injected", "depth_recovered", "mean_depth_final_lc", 'depth_relative_error',
                        "window_factor", "u1", "u2","r_star", "m_star", "loop_value"])
        print('The results are written in', self.final_filename)
        print("")
        return None
        
    def GetPlatoTimes(self):
        self.catalog_Q23 = []
        self.catalog_Q24 = []
        self.lcs_Q23 = []
        self.lcs_Q24 = []
        with open(f'../{self.ntbk_str}data/per_camera_data/c01_Q23_times.pkl', 'rb') as f:
            self.times_23 = pickle.load(f) ## 1D array
        with open(f'../{self.ntbk_str}data/per_camera_data/c01_Q24_times.pkl', 'rb') as f:
            self.times_24 = pickle.load(f)
        self.times_23 = self.times_23.reshape(self.initial_nb_points, self.bin_factor)
        self.times_23 = (np.mean(self.times_23, axis=1) - constants.plato_time_origin) / constants.seconds_per_day
        self.times_24 = self.times_24.reshape(self.initial_nb_points, self.bin_factor)
        self.times_24 = (np.mean(self.times_24, axis=1) - constants.plato_time_origin) / constants.seconds_per_day
        self.times = np.concatenate([self.times_23, self.times_24])
        self.void_lc = lk.LightCurve(time = self.times) 
        return None
        
    def GetPlatoFluxes(self):  
        for loop, camera in enumerate(prm.plato_cameras):
            if loop % 2 == 0:
                print("Camera", camera)
            camera_str = str(camera).zfill(2)
            if self.bin_factor == 1:
                self.file_name = f'{self.ntbk_str}data/per_camera_data/c{camera_str}'
            else:
                self.file_name = f'{self.ntbk_str}data/per_camera_data/binned{self.bin_factor}_c{camera_str}'
            
            with open(f'../{self.file_name}_Q23_fluxes.pkl', 'rb') as f:
                self.lcs_Q23.append(pickle.load(f))           ## 2D numpy array (one lightcurve per star) 
            with open(f'../{self.file_name}_Q24_fluxes.pkl', 'rb') as f:
                self.lcs_Q24.append(pickle.load(f))
            with open(f'../{self.ntbk_str}data/per_camera_data/c{camera_str}_Q23_catalog.pkl', 'rb') as f:
                self.catalog_Q23.append(pickle.load(f))       ## 2D array with 4 rows: position (0,1), mag (2) and ID (3)
            with open(f'../{self.ntbk_str}data/per_camera_data/c{camera_str}_Q24_catalog.pkl', 'rb') as f:
                self.catalog_Q24.append(pickle.load(f))       ## Columns are star IDs  
        print("GetData done")
        return [self.lcs_Q23, self.lcs_Q24, self.catalog_Q23, self.catalog_Q24, self.times]
    
    def GetPlatoLightCurves(self, data):
        self.stars_to_load = constants.plato_sunlike_stars
        self.nb_stars_used = len(self.stars_to_load)
        self.all_lcs = []
        self.stars_with_data = []
        self.has_Q23 = []
        self.has_Q24 = []
        self.lcs_Q23, self.lcs_Q24, self.catalog_Q23, self.catalog_Q24, self.times = data
        self.varsource = pd.read_feather(f'{self.ntbk_str}data/table_varsource_P1.ftr')
        self.varsource = self.varsource.fillna(0)
        self.stars_without_transit = np.where(pd.isnull(self.varsource['P_day']))[0]
        self.magnitudes = []
        for star_indice in range(constants.plato_total_nb_stars):
            if star_indice in self.stars_to_load:
                has_Q23_current = []
                has_Q24_current = []
                if star_indice % 200 == 0:
                    print("Star indice", star_indice)
                nb_23 = 0
                nb_24 = 0
                average_lc23 = np.zeros(self.initial_nb_points)
                average_lc24 = np.zeros(self.initial_nb_points)
                for loop in range(len(prm.plato_cameras)):
                    current_lc23 = None
                    current_lc24 = None
                    index_23 = np.where(np.abs(self.catalog_Q23[loop][3]-(star_indice+1))<0.01)[0]
                    index_24 = np.where(np.abs(self.catalog_Q24[loop][3]-(star_indice+1))<0.01)[0]
                    if len(index_23) > 0:
                        current_lc23 = self.lcs_Q23[loop][index_23[0]]
                        has_Q23_current.append("Yes")
                        current_lc23 = current_lc23.reshape(self.initial_nb_points)
                        average_lc23 += current_lc23
                        nb_23 += 1
                        mag = self.catalog_Q23[loop][2, index_23[0]]
                    else: 
                        has_Q23_current.append("No")
                        current_lc23 = np.empty(1)
                        mag = 0
                    if len(index_24) > 0:
                        current_lc24 = self.lcs_Q24[loop][index_24[0]]
                        has_Q24_current.append("Yes")
                        current_lc24 = current_lc24.reshape(self.initial_nb_points)
                        average_lc24 += current_lc24
                        nb_24 += 1
                        mag = self.catalog_Q24[loop][2, index_24[0]]
                    else: 
                        has_Q24_current.append("No")
                        current_lc24 = np.empty(1)
                        mag = max(0, mag)
                if prm.is_notebook:
                    self.show_data = prm.initial_plot_data
                if nb_23 > 0 and nb_24 > 0:
                    average_lc = np.concatenate([average_lc23/nb_23, average_lc24/nb_24])
                    t = self.times
                elif nb_23 > 0:
                    average_lc = average_lc23 / nb_23
                    t = self.times_23
                elif nb_24 > 0:
                    average_lc = average_lc24 / nb_24
                    t = self.times_24
                else: 
                    average_lc = self.void_lc 
                    if prm.is_notebook:
                        self.show_data = False
                if prm.is_notebook:
                    if self.show_data:
                        lc = lk.LightCurve(time=t, flux=average_lc);
                        lc.plot();
            else:
                mag, has_Q23_current, has_Q24_current, average_lc = None, None, None, None

            self.magnitudes.append(mag)
            self.has_Q23.append(has_Q23_current)
            self.has_Q24.append(has_Q24_current)  
            self.all_lcs.append(average_lc)

        return [self.has_Q23, self.has_Q24, self.all_lcs, self.stars_without_transit]
    
    def GetPlatoDataQuick(self, all_data):
        self.times, self.has_Q23, self.has_Q24, self.all_lcs, self.stars_without_transit = all_data
        self.stars_without_transit = self.stars_without_transit.tolist()
        self.varsource = pd.read_feather(f'../{self.ntbk_str}data/table_varsource_P1.ftr').fillna(0)
        self.DumpSunLike()
        return None

    def DumpSunLike(self):
        name = f'../{self.ntbk_str}data/all-sunlike_both-lcs_bin{self.bin_factor}.pkl'
        print("File name:", name)
        with open (name, 'wb') as x:
                pickle.dump([self.magnitudes, self.has_Q23, self.has_Q24, self.all_lcs, self.stars_without_transit], x)
        return None

    def GetPlatoSunlike(self):
        name = f'../{self.ntbk_str}data/all-sunlike_both-lcs_bin{self.bin_factor}.pkl'
        with open(name, 'rb') as x:
            self.magnitudes, self.has_Q23, self.has_Q24, self.all_lcs, self.stars_without_transit = pickle.load(x)
        self.varsource = pd.read_feather(f'../{self.ntbk_str}data/table_varsource_P1.ftr').fillna(0)
        return None

    def GetKeplerDataQuick(self):
        self.star_ids = prm.kepler_stars
        self.nb_stars_used = len(self.star_ids)
        with open(f'../{self.ntbk_str}data/kepler_lightcurves_Q2-Q3.pkl', 'rb') as file:
            self.times, self.all_lcs, self.timestamps = pickle.load(file)
        self.nb_data_points = len(self.times)
        return None

    def DurationGrid(self):
        self.bls_min_dur = max(prm.bls_duration_extremes[0], 0.7*(prm.bls_minimum_period/constants.days_in_year)**(1/3)*constants.earth_transit_dur)
        self.bls_max_dur = min(prm.bls_duration_extremes[1], 1.2*(prm.bls_maximum_period/constants.days_in_year)**(1/3)*constants.earth_transit_dur)
        if prm.transit_search == "bls":
            self.bls_duration_grid = []
            dur = self.bls_min_dur
            while dur < self.bls_max_dur:
                self.bls_duration_grid = self.bls_duration_grid + [dur]
                dur *= prm.bls_duration_grid_step
            if dur / self.bls_max_dur < prm.bls_duration_grid_step / 1.2:
                self.bls_duration_grid = self.duration_grid + [dur]
        return None
        
    def BeforeLoop(self):
        print("Mission:", prm.mission)
        print("Transit search:", prm.transit_search)
        if prm.mission == 'plato':
            print("Cadence:", prm.plato_cadence, "seconds")
            if prm.plato_stars == constants.plato_sunlike_stars:
                print("Using all PLATO 'sun-like stars', and no other")
            if self.nb_stars_used == 1:
                self.star_mag, self.star_nb_cams = constants.star_params[prm.plato_stars[0]][2:4]
                print("Star apparent magnitude:", self.star_mag)
        if prm.is_notebook:
            prm.nb_samples = 1
            if prm.reproduce_csv==False:
                print('Using random values')
            else: 
                print('Reproducing results')
                prm.nb_samples = 1  
        print("")
        print("Numpy seed:", prm.numpy_seed)
        print("Detrending before injection?", prm.detrend_before_injection)
        print("Planet radius bounds:", prm.min_radius, prm.max_radius)
        if prm.mission == 'plato':
            self.GetPlatoTimes()
            self.GetPlatoSunlike()
            if prm.create_new_lc_files:
                data = self.GetPlatoFluxes()
                all_data = [data[-1]] + self.GetPlatoLightCurves(data)
                self.GetPlatoDataQuick(all_data)
        else:
            self.GetKeplerDataQuick()
        if prm.is_notebook == False:
            self.OpenWriter()
            if prm.nb_samples > 49:
                shutil.copyfile('current_parameters.py', f'../results/{self.final_filename}_params.py')
        self.loop_value = 0
        self.count_samples = 0
        self.previous_time = time.time()
        self.DurationGrid()
        return None


####################################### PART 2: functions inside of the loop


    def RandomLaws(self, law, min_value, max_value):
        if law == "uniform":
            return np.random.uniform(min_value, max_value)
        elif law == "log-uniform":
            return min_value * (max_value / min_value) ** np.random.random()
        elif law == "constant":
            return min_value

    def RandomValues(self): ##Changed order
        # Accepted
        self.planet_radius = constants.earth_sun_ratio * self.RandomLaws(law=prm.radius_law, min_value=prm.min_radius, max_value=prm.max_radius)
        if prm.transit_search == "gpu":
            self.true_period = self.RandomLaws(law=prm.monos_period_law, min_value=prm.monos_min_period, max_value=prm.monos_max_period)
        else:
            self.true_period = self.RandomLaws(law=prm.multis_period_law, min_value=prm.multis_min_period, max_value=prm.multis_max_period)
        self.impact_parameter = self.RandomLaws(law=prm.impact_law, min_value=prm.min_impact, max_value=prm.max_impact)
        self.true_t0 = self.true_period*np.random.random()

        self.i = int(np.random.random()*self.nb_stars_used) 
        self.m_star, self.r_star = 1, 1
        if prm.mission == "plato":
            self.star_key = prm.plato_stars[self.i]
            self.star_index = constants.star_params[self.star_key][0]
            self.star_mag, self.star_nb_cams = constants.star_params[self.star_key][2:4]
            self.chosen_quarter = constants.star_params[self.star_key][1]
            if prm.true_star_params:
                self.m_star = self.varsource.iloc[self.star_index]['Ms_Msun']
                self.r_star = self.varsource.iloc[self.star_index]['Rs_Rsun']
        else:
            self.star_index = self.star_ids[self.i]
        return None

    def ReproduceValues(self):
        df = pd.read_csv(prm.csv_file)
        self.true_period = df.loc[self.n,'period_days']
        self.planet_radius = df.loc[self.n,'Radius (r_earth)']
        self.impact_parameter = df.loc[self.n,'Impact parameter']
        self.true_t0 = df.loc[prm.csv_row_number,'t0_days']
        self.star_index = df.loc[prm.csv_row_number,'Star ID']
        self.m_star = self.varsource.iloc[self.star_index]['Ms_Msun']
        self.r_star = self.varsource.iloc[self.star_index]['Rs_Rsun']
        print("R_star ok")
        return None

    def TransitType(self):
        if self.impact_parameter < 1 - self.planet_radius/self.r_star:
            self.transit_type = "Full transit"
        elif self.impact_parameter < 1 + self.planet_radius/self.r_star:
            self.transit_type = "Partial transit"
        else:
            self.transit_type = "No transit"
        return None
    
    def PlatoTimeArray(self, period):
        n_tot = np.floor((constants.plato_Q24_end-self.true_t0)/period)+1
        n0 = np.floor((constants.plato_Q24_start-self.true_t0)/period)+1
        n1 = np.floor((constants.plato_Q23_end-self.true_t0)/period)+1
        n2 = n_tot - n0
        self.quarter_message = ""
        if "Yes" in self.has_Q23[self.star_index] and ("Yes" in self.has_Q24[self.star_index]) and self.chosen_quarter=="both":
            nb_transits = n_tot
            self.t = self.times
            self.allow_cut = True
        elif "Yes" in self.has_Q23[self.star_index]:
            nb_transits = n1
            self.t = self.times_23
            self.allow_cut = False
            if self.chosen_quarter == "both":
                self.quarter_message = "Careful, we are only using Q23"
            elif self.chosen_quarter == "Q24":
                raise ValueError("You asked for Q24 but the star only has Q23")
        elif "Yes" in self.has_Q24[self.star_index]:
            nb_transits = n2
            self.t = self.times_24
            self.allow_cut = False
            if self.chosen_quarter == "both":
                self.quarter_message = "Careful, we are only using Q24"
            elif self.chosen_quarter == "Q23":
                raise ValueError("You asked for Q23 but the stars only has Q24")
        else: 
            self.transit_type = "No data"
            nb_transits = None
        return nb_transits

    def ComputeKeplerTransits(self, period):
        n_tot = np.floor((constants.kepler_Q3_end-self.true_t0)/self.true_period)+1
        n0 = np.floor((constants.kepler_Q2_start-self.true_t0)/self.true_period)+1
        nb_transits = n_tot - n0
        return nb_transits

    def BeforeModel(self):
        self.loop_value += 1
        if prm.is_notebook and prm.reproduce_csv:
            self.ReproduceValues()
        else:
            self.RandomValues()
        self.TransitType()
        if prm.mission == "plato":
            self.nb_transits_injected = self.PlatoTimeArray(period=self.true_period)
        else:
            self.nb_transits_injected = self.ComputeKeplerTransits(period=self.true_period)
        return None

    def CarryLoop(self):
        return self.loop_value >= prm.min_loop_value and (self.transit_type == "Full transit" or (self.transit_type == "Partial transits" and prm.include_partial_transits)) and self.nb_transits_injected >= self.min_nb_transits and self.nb_transits_injected <= self.max_nb_transits

    def CreateModel(self):
        self.count_samples += 1
        if prm.mission == "kepler":
            self.t = np.array(self.times[0])
        params = batman.TransitParams()
        params.t0 = self.true_t0                   #time of inferior conjunction
        params.per = self.true_period                 #orbital period in days
        params.rp = self.planet_radius/self.r_star                  #planet radius (in units of stellar radii)
        params.a = ((params.per*constants.seconds_per_day)**2 * constants.G * (self.m_star*constants.M_sun) / (4*np.pi**2))**(1/3) / (constants.R_sun*self.r_star) #semi-major axis (in units of stellar radii)
        params.inc = 180/np.pi*np.arccos(self.impact_parameter/params.a)          #orbital inclination (in degrees)
        params.ecc = 0.                    #eccentricity
        params.w = 90.                       #longitude of periastron (in degrees)
        params.u = prm.u                #limb darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"       #limb darkening model
        self.a = params.a
        self.model_batman = batman.TransitModel(params, self.t)
        self.array_model = self.model_batman.light_curve(params)
        self.lc_model = lk.LightCurve(time=self.t, flux=self.array_model)
        if prm.is_notebook: 
            if prm.plot_model:
                self.lc_model.plot();
        return None

    def NoTransit(self):
        if prm.mission == "kepler":
            self.array_data = self.all_lcs[self.i]
        else:
            self.array_data = self.all_lcs[self.star_index]
        self.lc_result = lk.LightCurve(time=self.t, flux=self.array_data)
        return None

    def JustTransitModel(self):
        self.lc_result = lk.LightCurve(time=self.t, flux=self.array_model)

    def ReducedNoise(self):
        if prm.mission == "kepler":
            self.array_data = self.all_lcs[self.i]
        else:
            self.array_data = self.all_lcs[self.star_index]
        
        self.lc_noise = lk.LightCurve(time=self.t, flux=self.array_data)
        self.window_length = 2*int(prm.window_factor/2 * self.bls_max_dur/prm.plato_Q23_end*self.initial_nb_points) + 1

        if self.count_samples == 1:
            self.flat_noise = self.lc_noise.flatten(window_length = self.window_length)
            self.ShowDuration(message = "pre-injection flattening duration:")
        
        self.array_result = (prm.noise_reduction*(self.flat_noise.flux-1)+1) * self.array_model
        self.lc_result = lk.LightCurve(time=self.t, flux=self.array_result)
        return None
            
    def TransitInjection(self):
        if prm.mission == "kepler":
            self.array_data = self.all_lcs[self.i]
        else:
            self.array_data = self.all_lcs[self.star_index]
        self.array_result = self.array_data * self.array_model
        self.lc_result = lk.LightCurve(time=self.t, flux=self.array_result)
        return None

    def Injection(self):
        if prm.include_transit == False:
            self.NoTransit()
        elif prm.include_data == False:
            self.JustTransitModel()
        elif prm.detrend_before_injection:
            self.ReducedNoise()
        else:
            self.TransitInjection()
        return None

    def FlatteningRemoval(self):
        if prm.mission == 'plato':
            if self.chosen_quarter == "Q23" and self.allow_cut:
                self.lc_result = self.lc_result[:len(self.lc_result.flux)//2]
            elif self.chosen_quarter == "Q24" and self.allow_cut:
                self.lc_result = self.lc_result[len(self.lc_result.flux)//2:]
            self.window_length = 2*int(prm.window_factor/2 * self.bls_max_dur/constants.plato_Q23_end*self.initial_nb_points) + 1
        else:
            self.initial_nb_points = len(self.t)
            self.window_length = 2*int(prm.kepler_window_factor/2 * self.bls_max_dur/constants.plato_Q23_end*self.initial_nb_points) + 1
        
        if prm.detrend_before_injection:
            if prm.remove_outliers:
                self.final_lc = self.lc_result.remove_outliers(sigma_upper=prm.sigma_upper, sigma_lower=prm.sigma_lower*self.bin_factor**0.5)
            else:
                self.final_lc = self.lc_result
            return None
        self.flat_lc = self.lc_result.flatten(window_length = self.window_length)
        if prm.remove_outliers:
            self.final_lc = self.flat_lc.remove_outliers(sigma_upper=prm.sigma_upper, sigma_lower=prm.sigma_lower)
        else:
            self.final_lc = self.flat_lc
        self.flat_duration = self.ShowDuration()
        self.sum_computing_detrending += self.flat_duration
        self.t = self.final_lc.time.value
        self.final_nb_points = len(self.t)
        self.array_fluxes = self.final_lc.flux.value
        self.star_noise = 10**6*np.std(self.array_fluxes)
        return None
    
    def ShowCurves(self, which_ones):
        self.title = "Sample n°"+str(self.count_samples)
        self.xlabel = "Sample n°"+str(self.count_samples)+" Period is "+str(self.true_period)
        self.ylabel = "Radius^2 is "+str(np.round(self.planet_radius**2,5))
        if 0 in which_ones:
            self.lc_result.plot(title=self.title, xlabel=self.xlabel, ylabel=self.ylabel);
        if 1 in which_ones:
            self.flat_lc.plot(title=self.title, xlabel=self.xlabel, ylabel=self.ylabel);
        if 2 in which_ones:
            self.final_lc.plot(title=self.title, xlabel=self.xlabel, ylabel=self.ylabel);
        return None

    def BlsSearch(self):
        self.positive_power_threshold = prm.bls_power_threshold / self.bin_factor ** 0.5
        self.SNR_threshold = prm.bls_SNR_threshold
        self.bls = self.final_lc.to_periodogram('bls', frequency_factor=prm.bls_frequency_factor, 
                                                minimum_period=prm.bls_minimum_period, duration=self.bls_duration_grid,
                                               maximum_period=prm.bls_maximum_period)   
        self.measured_period = self.bls.period_at_max_power.value  
        self.period_relative_error = -(self.true_period-self.measured_period)/self.true_period
        self.measured_t0 = self.bls.transit_time_at_max_power.value
        self.measured_dur = self.bls.duration_at_max_power.value
        self.measured_depth = self.bls.depth_at_max_power.value*10**6
        self.search_power = self.bls.max_power
        
        self.bls_duration = self.ShowDuration()
        self.sum_computing_transit_search += self.bls_duration
        if prm.is_notebook:
            if prm.plot_periodogram:
                self.bls.plot();
        return None

    def TlsSearch(self):
        self.positive_power_threshold = prm.tls_positive_power_threshold
        self.SNR_threshold = prm.bls_SNR_threshold
        tls_model = transitleastsquares(self.t, self.array_fluxes)
        self.tls_results = tls_model.power(oversampling_factor=prm.tls_oversampling_factor, duration_grid_step = prm.tls_duration_grid_step)
        self.measured_period = self.tls_results.period
        self.period_relative_error = -(self.true_period-self.measured_period)/self.true_period
        self.measured_t0 = self.tls_results.T0
        self.measured_dur = self.tls_results.duration
        self.measured_depth = (1-self.tls_results.depth)*10**6
        self.search_power = self.tls_results.SDE
        self.bls_duration = self.ShowDuration()
        self.sum_computing_transit_search += self.bls_duration
        return None

    def GPUTransitSearch(self):
        self.positive_power_threshold = prm.gpu_power_threshold
        self.SNR_threshold = prm.gpu_SNR_threshold
        # add errors
        if prm.include_data==False:
            self.errors_value = 0.0001
        else:
            self.errors_value = prm.gpu_error_factor * np.median(np.abs(np.median(self.array_fluxes) - self.array_fluxes))
        self.errors = np.repeat(self.errors_value, self.final_nb_points)

        # transit duration range to check
        min_duration = 0.04  # ~1 hour
        max_duration = 1.00  # 24 hours
        
        # do the jit compilation for a fair comparison
        tced = TCEDetector(self.t, self.array_fluxes, self.errors,
            min_duration, max_duration, bin_width=None)
        
        # find the transits
        self.transits, chisq, self.depth = tced.find_transits(depth_threshold=prm.gpu_depth_threshold)
        
        self.measured_period = self.true_period
        self.period_relative_error = 0
        self.bls_duration = self.ShowDuration()
        self.sum_computing_transit_search += self.bls_duration
        
        keys = ['t_start_lsq', 'peak_gridsearch_power', 'duration_lsq','depth_lsq_ppm', 't_start_gridsearch', 'duration_gridsearch', 'depth_gridsearch']
        max_power = 0
        for transit in self.transits:
            new_power = transit[keys[1]]
            if new_power > max_power:
                powerful_transit = transit
                max_power = new_power
        if keys[2] in powerful_transit:
            self.measured_dur = powerful_transit[keys[2]]
            self.measured_t0 = powerful_transit[keys[0]] + self.measured_dur/2
            self.measured_depth = powerful_transit[keys[3]]
        else: 
            self.measured_dur = powerful_transit[keys[5]]
            self.measured_t0 = powerful_transit[keys[4]]
            self.measured_depth = powerful_transit[keys[6]]
        self.search_power = max_power
        return None

    def TransitSearch(self):    
        if prm.transit_search == 'bls':
            self.BlsSearch()
        elif prm.transit_search == 'tls':
            self.TlsSearch()
        elif prm.transit_search == 'gpu':
            self.GPUTransitSearch()
        else:
            raise ValueError("Transit search ill-defined: choose 'bls', 'tls' or 'gpu'")
        return None

    def MeanDepth(self):
        measured_mean_depth = 0
        start_time = (self.measured_t0 - self.measured_dur/2)
        end_time = start_time + self.measured_dur
        x = 2
        if prm.mission == 'plato':
            quarter_in_days = constants.plato_Q23_end
            if self.chosen_quarter != "both":
                x = 1
            if self.chosen_quarter == "Q24":
                start_time, end_time = start_time - constants.plato_Q24_start, end_time - constants.plato_Q24_start
        else:
            quarter_in_days = constants.kepler_Q2_end - constants.kepler_Q2_start
        self.nb_transits_recovered = 0
        add = True
        while start_time < self.t[-1] and (self.nb_transits_recovered < 1 or prm.transit_search != "gpu"):
            if start_time > 89 and add:
                start_time -= 2
                end_time -= 2
                add = False
            start_point, end_point = int(start_time*self.final_nb_points/quarter_in_days/x), int(end_time*self.final_nb_points/quarter_in_days/x)+1
            measured_mean_depth += (1 - np.mean(self.final_lc.flux.value[start_point:end_point+1])) * 10**6
            start_time += self.measured_period
            end_time += self.measured_period
            self.nb_transits_recovered += 1
        self.measured_mean_depth = measured_mean_depth / self.nb_transits_recovered
        self.SNR = self.measured_depth/self.star_noise*(self.nb_transits_recovered*self.final_nb_points/quarter_in_days*self.measured_dur)**0.5
        return None

    def ShowDuration(self, save=False, message=None):
        self.current_time = time.time()
        computing_duration = self.current_time - self.previous_time
        if message != None:
            print(message, computing_duration)
        self.previous_time = self.current_time
        return computing_duration
    
    def FoldedCurve(self, whole_period):
        max_phase = prm.initial_folding_max_phase
        print("Bls power:", self.bls.max_power)
        if prm.folding_parameters == "measured":
            period = self.measured_period
            t0 = self.measured_t0
            dur = self.measured_dur
            print("Using measured values for folding")
        elif prm.folding_parameters == "true": 
            period = self.true_period
            t0 = self.true_t0
            dur = self.duration_grid[2]
        else: raise Exception("Please chose between 'measured' and 'true' for prm.folding_parameters")
        self.planet_b_model = self.bls.get_transit_model(period=period,
                                           transit_time=t0,
                                           duration=dur)
        ax = self.final_lc.fold(period, self.measured_t0).scatter()    
        self.planet_b_model.fold(period, self.measured_t0).plot(ax=ax, c='r', lw=1)    
        if whole_period: max_phase = self.measured_period/2
        ax.set_xlim(-max_phase,max_phase);
        return None
    
    def ComputeResults(self):
        self.t0_error = np.abs((self.measured_t0 - self.true_t0)%self.true_period)
        self.t0_error = min(self.t0_error, np.abs(self.true_period-self.t0_error))
        
        x = (1+self.planet_radius/self.r_star)**2 - self.impact_parameter**2
        if x>0:
            x1 = x**0.5/self.a
            if x1 > 1:
                print("Arcsin value problem:", x1)
            self.tdur = np.arcsin(x1)*self.true_period/np.pi

        # Limb darkening, and resulting injected transit depth
        mu = (1-self.impact_parameter**2)**0.5
        peak_darkening = 1 - prm.u[0]*(1-mu) - prm.u[1]*(1-mu)**2
        total_darkening = 1 - prm.u[0]/3 - prm.u[1]/6
        b = self.impact_parameter
        chord_darkening = 1 - prm.u[0]*(1 - np.pi/4*(1-b**2)**0.5) - prm.u[1]*(5/3 - np.pi/2/(1-b**2)**0.5-b**2*(2/3-np.pi/2/(1-b**2)**0.5))
        self.theo_max_depth = (self.planet_radius/self.r_star)**2*peak_darkening*10**6/total_darkening
        self.theo_avg_depth = (self.planet_radius/self.r_star)**2*chord_darkening*10**6/total_darkening
        
        # Recovered depth and dur
        if prm.transit_search == "gpu":
            self.depth_relative_error = self.measured_depth/self.theo_max_depth
        else:
            self.depth_relative_error = self.measured_depth/self.theo_avg_depth
        self.dur_relative_error = self.measured_dur/self.tdur

    def Outcomes(self):
        harm = 0
        is_detectable = self.nb_transits_injected > 0 and (self.nb_transits_injected > 1 or prm.transit_search == 'gpu')
        if self.nb_transits_injected > 0:
            self.count_transits += 1
        # Different outcomes
        if self.nb_transits_injected > 0 and (np.abs(self.period_relative_error) < prm.max_relative_period_error or (self.measured_t0 - self.measured_period < self.t[0] and self.measured_t0 + self.measured_period > self.t[-1])) and self.t0_error < prm.max_t0_error:
            if self.SNR > self.SNR_threshold and self.search_power > self.positive_power_threshold:
                harm = 100
                self.correct_average += 1
                if self.nb_transits_injected == 1 or prm.transit_search == 'gpu':
                    self.result_message = "True mono-transit positive"
                    self.outcome = "true_pos_mono"
                else:
                    self.result_message = "True multi-transits positive"
                    self.outcome = "true_pos_multi"
                self.depth_ratios = self.depth_ratios + [self.measured_depth/self.theo_max_depth]
                self.avg_depth_ratios = self.avg_depth_ratios + [self.measured_depth/self.theo_avg_depth]
                self.dur_ratios = self.dur_ratios + [self.dur_relative_error]
            else:
                self.result_message = "FALSE NEGATIVE RELATED TO SIGNAL!"
                self.outcome = "false_neg_related"
        elif self.nb_transits_injected > 0 and (np.abs(self.period_relative_error + 0.5) < 0.01 or np.abs(self.period_relative_error - 1) < prm.max_relative_period_error or np.abs(self.period_relative_error - 2) < prm.max_relative_period_error) and self.t0_error < prm.max_t0_error: 
            if self.SNR > self.SNR_threshold and self.search_power > self.positive_power_threshold:
                harm = 100   
                self.result_message = "HARMONIC POSITIVE"
                self.outcome = "harmonic_pos"
            else:
                self.result_message = "FALSE NEGATIVE RELATED TO SIGNAL (harmonic)!"
                self.outcome = "false_neg_related"
                self.below_thresholds_avg += 1
        else:
            if self.SNR > self.SNR_threshold and self.search_power > self.positive_power_threshold: 
                self.result_message = "FALSE POSITIVE!"
                self.outcome = "false_pos"
                self.false_pos_average += 1
            elif self.nb_transits_injected > 0:
                self.result_message = "False negative, unrelated to signal"
                self.outcome = "false_neg_unrelated"
            else:
                self.result_message = "True negative: no transit injected"
                self.outcome = "true_pos"
        return None

    def WriteResults(self):
        self.writer.writerow([prm.numpy_seed, self.count_samples, self.star_index, self.chosen_quarter, self.star_mag, self.star_nb_cams,
                              prm.transit_search, self.planet_radius, self.impact_parameter, self.transit_type,
                              self.nb_transits_injected, self.nb_transits_recovered, self.true_period, self.measured_period, self.period_relative_error, 
                              self.true_t0, self.measured_t0, self.t0_error, self.correct_average*100/max(1,self.count_transits), 
                              self.star_noise, self.final_nb_points, self.SNR, self.search_power, self.outcome, 
                              self.tdur, self.measured_dur, self.dur_relative_error,
                              self.theo_max_depth, self.theo_avg_depth, self.measured_depth, self.measured_mean_depth, self.depth_relative_error,
                              prm.window_factor, prm.u[0], prm.u[1], self.r_star, self.m_star, self.loop_value])   
        return None

    def LoopOutputs(self):
        ## Main infos
        print("Computing sample n°",self.count_samples, "of", prm.nb_samples)
        print(self.result_message)
        print("")
        print("Transit search:", prm.transit_search)
        if prm.detrend_before_injection:
            print("Detrending before injection!")
        if prm.remove_outliers == False:
            print("Remove outliers?", prm.remove_outliers)
        if self.transit_type != "Full transit":
            print(self.transit_type)
        print('Star_index', self.star_index)
        print("Quarter:", self.chosen_quarter)
        print("Number of transits injected:", int(self.nb_transits_injected))
        print("Number of transits recovered:", self.nb_transits_recovered)
        if prm.mission == "plato" and self.bin_factor == 1:
            print("Using unbinned Plato curves")
        print("")

        ## Period and t0
        print("True period:", self.true_period)
        if prm.transit_search != 'gpu':
            print("Measured period:", self.measured_period)
            print("Relative error:", self.period_relative_error)
        if self.quarter_message != "":
            print(self.quarter_message)
        print("True t0:", self.true_t0)
        print("Measured t0:", self.measured_t0)
        print("t0 error:", self.t0_error)
        print("")

        ## Depth
        print("Planet radius in r_earths", self.planet_radius/constants.earth_sun_ratio)
        print("Impact parameter is", self.impact_parameter)
        print('Theoretical max depth in ppm', self.theo_max_depth)
        print('Theoretical mean depth in ppm', self.theo_avg_depth)
        print('Measured depth in ppm', self.measured_depth)
        print('Measured mean depth', self.measured_mean_depth)
        if self.measured_depth != None:
            if self.nb_transits_injected == 1:
                print("Depth ratio", self.measured_depth/self.theo_max_depth)
            else:
                print("Depth ratio", self.measured_depth/self.theo_avg_depth)
        print("")
        print("Search power:", self.search_power)
        print("SNR:", self.SNR)
        print('Star std', self.star_noise)
       
        ## Transit duration
        print("Theoretical duration", self.tdur)
        print("Measured duration", self.measured_dur)
        print("Duration ratio is", self.measured_dur/self.tdur)
        print("")

        ##ShowDuration
        if prm.transit_search == 'Bls' and self.nb_transits_injected > 1:
            print("Frequency factor is", prm.bls_frequency_factor)
        print("Transit search duration", self.bls_duration)
        print("Detrending duration", self.flat_duration)
        print("Average transit search duration :", self.sum_computing_transit_search/max(1, self.count_samples))
        print("Average detrending duration :", self.sum_computing_detrending/max(1, self.count_samples))
        print("")
        
        if prm.print_all_values:
            print("Darkening coefficients", prm.u)
            print('Sigma lower and upper are', prm.sigma_lower, prm.sigma_upper)
            if prm.transit_search == 'bls' and self.nb_transits_injected > 1:
                print('Min and max bls period:', prm.bls_minimum_period, prm.bls_maximum_period)
                print('Speed factor', self.speed_factor)
                print("Frequency factor is", prm.bls_frequency_factor)
            print("Star mass is", self.m_star)
            print("Star radius is", self.r_star)
            if prm.mission == "plato":
                print("Using cameras", prm.plato_cameras)

        ##Flatremove
        if prm.print_all_values:
            print('Window length is', self.window_length, "or", self.window_length/self.initial_nb_points*constants.plato_Q23_end, "days") 

        ##Writerseults
        print("Current true positive rate is", self.correct_average*100/max(1,self.count_transits))
        print("Current false positive rate is", self.false_pos_average*100/max(1,self.count_samples))
        print("Current false negative linked to injected signal (below thresholds) rate:", self.below_thresholds_avg*100/max(1,self.count_transits))
        if prm.is_notebook == False:
            print('The results were written in', self.final_filename)
        print("")
        print("")
        return None

    def AfterLoop(self):
        print("Final computing averages transit search:", self.sum_computing_transit_search/max(self.count_samples,1))
        print("Final computing averages detrending:", self.sum_computing_detrending/max(self.count_samples,1))
        print("")
        if len(self.depth_ratios) > 0:
            print('Mean and std of max depth ratios (should be 1 for gpu)', np.mean(self.depth_ratios), np.std(self.depth_ratios))
            print('Mean and std of mean depth ratios (should be 1 for bls)', np.mean(self.avg_depth_ratios), np.std(self.avg_depth_ratios))
            print('Mean and std of dur ratio', np.mean(self.dur_ratios), np.std(self.dur_ratios))
        if prm.is_notebook == False:
            self.fgh.close()
        print("")
        print("")
        print("")
        return None

