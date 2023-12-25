import constants

######################## main
# Any positive integer; ensures that results can be reproduced:
numpy_seed = 0

# Number of transits to be injected and searched:
nb_samples = 10

# To start at a specific loop value:
min_loop_value = 0

# If true, more details are shown when running program
print_all_values = False



######################## get_data
# "plato" or "kepler"
mission = "plato"                  

### If "plato"

# 25, 100 or 300 (seconds)
plato_cadence = 100

# Format: "{star_id}_{quarter}"
plato_stars = ["624_both", "790_both"] 

# If needed, only include certain cameras; format: list of integers, between 1 and 24
plato_cameras = constants.all_plato_cameras

## If "kepler"; format: list of integers
kepler_mag = 10
kepler_stars = constants.kepler_sunlike[kepler_mag]                           



######################## model_transit
# Choose radius (in r_earths), period (in days) and impact parameter bounds
min_radius, max_radius = 0.5, 4
min_impact, max_impact = 0, 0.2
multis_min_period, multis_max_period = 45, 180
monos_min_period, monos_max_period = 45, 180

# Choose "constant", "uniform" or "log-uniform"; 
# If "constant", all models use min_radius or min_period value
radius_law = "log-uniform"
impact_law = "uniform"
multis_period_law = "log-uniform"
monos_period_law = "log-uniform"

# If false, only monos for gpu, and all nb_transits >= 2 for bls and tls
override_nb_transits = True
# Only relevant if override is true; bounds included
min_nb_transits, max_nb_transits = 0, constants.inf                                        
# If false, samples are only modeled if impact_parameter + r_planet / r_star < 1 (the planet comes fully in front of the star)
include_partial_transits = False

# Quadratic limb darkening coefficients
u = [0.48, 0.19]
# If false, all 'sun-like stars' are considered to have R=R_sun and M=M_sun (this makes the results easier to analyze)
true_star_params = False  



######################## inject_detrend
# If false, transit search is performed directly on the transit model, to check that it works ok
include_data = True
# If false, the search is performed on the Plato/Kepler data with no injected transit -> can be interesting to test false positives
include_transit = True

# If true, the Plato/Kepler lightcurve is detrended before transit injection; unrealistic but useful for comparison
# In that case, 'noise reduction' allows to reduce (if set to <1 value) the Kepler/Plato noise, also intersting for comparison
detrend_before_injection = False
noise_reduction = 1

remove_outliers = True
sigma_lower = 6
sigma_upper = 5
window_factor = 15
kepler_window_factor = 5

noise_reduction = 1 
raw_curves = False                                                 



######################## transit_search
transit_search = "bls"

# If "tls"
tls_oversampling_factor = 1
tls_duration_grid_step = 1.25

# If "bls"
# For 100s bins, min_period=1.2, max_period=100, frequency_factor is good at 20-25
bls_frequency_factor = 50
bls_minimum_period = 1.5
bls_maximum_period = 185

bls_duration_grid_step = 1.2
bls_duration_extremes = [0.035, 1]
bls_adapt_speed = False

# If "gpu"
gpu_depth_threshold = None
gpu_error_factor = 1.4826



#################################### Results
# Two files are created: file_name.csv, containing all results, and params_file_name.py, which is a copy of this file with the parameters used
# If 'default', the file names will be chosen to display the main informations,
# and to ensure they are new (so that the previous files are not erased)
base_filename = 'default'

max_relative_period_error = 0.01
max_t0_error = 0.25

bls_SNR_threshold = 7
# Defined for 25s cadence, and scaled by a factor 1 / sqrt(cadence) for higher ones
bls_power_threshold = 1  * 10**(-5)
bls_kepler_power_threshiold = 1 * 10**(-6)

tls_positive_power_threshold = 9

gpu_SNR_threshold = 6
gpu_power_threshold = 14


#################################### Others
reproduce_csv = False
is_notebook = False
create_new_lc_files = False

# Plots
folding_parameters = "measured"
initial_folding_max_phase = 1
initial_plot_data = False
plot_model = False
plot_periodogram = False                       

# Reproduce
csv_file = f'../csvs/plato_csvs/bin4_seed0.csv'
csv_row_number = 956 
notebook_string = "../exoplanets_internship/"
