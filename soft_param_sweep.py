import numpy as np
import mwahpy.output_handler
import matplotlib.pyplot as plt
import os.path
import subprocess
import pandas as pd

def set_vars(mass, rscale, model, r1, rc):
    vol_pcrit = 0.568910904587397184785763397846734505212216314432372653620
    pcrit = 0.000679087369829744220469326744094105320596648627735869652
    r200 = (mass / vol_pcrit) ** (1.0 / 3.0)

    if model == 'cored':
        D1 = r1 * (1 + r1 / rscale)**2 / (rscale + rscale * (r1 / rc)**2)
        D2 = (rscale)**3 * (np.log(1 + r200 / rscale) - np.log(1 + r1 / rscale) - r200 / (rscale + r200) + r1 / (rscale + r1))
        D3 = (rc)**2 * (r1 / (1 + (rc / r1)**2) - rc * np.arctan(r1 / rc) + r1 / (1 + (r1 / rc)**2))
        p0 = mass / (4 * np.pi * (D1 * D2 + D3))
        ps = p0 * D1
    elif model == "nfw":  # NFW
        c = r200 / rscale
        term = np.log(1.0 + c) - (c / (1.0 + c))
        p0 = 200.0 * c**3 * pcrit / (3.0 * term)
        ps = None

    bound = 5.0 * r200

    return r200, p0, ps, bound

def get_extra_nfw_mass(p0, ps, bound, model, rscale, r1, rc):
    r = bound
    rs = rscale
    if model == 'cored':
        C1 = 0
        C3 = C1 + 4 * np.pi * (p0 * (rc**2) * ((r1**3) / ((r1**2) + (rc**2)) - rc * np.arctan(r1 / rc) + r1 / (1 + (r1 / rc)**2)) - ps * (rs**3) * (np.log(1 + r1 / rs) - r1 / (rs + r1)))
        if r <= r1:
            m = 4.0 * np.pi * p0 * (rc**2) * (r / (1 + (rc / r)**2) - rc * np.arctan(r / rc) + r / (1 + (r / rc)**2)) - C1
        else:
            m = 4.0 * np.pi * (rs**3) * ps * (np.log(1 + r / rs) - r / (rs + r)) - C3
    elif model == "nfw":  # NFW
        m = 4.0 * np.pi * p0 * (rs**3) * (np.log((rs + r) / rs) - r / (rs + r))
    return m

def density_function(scale_length, mass, radius, model, r1=None, rc=None, p0=None, ps=None):
    if model == "plummer":
        density = (
            (3 * mass)
            / (4 * np.pi * scale_length**3)
            * (1 + radius**2 / scale_length**2) ** (-5 / 2)
        )
    elif model == "cored":
        rs = scale_length
        density = []
        for r in radius:
            if r <= r1:
                density.append( 
                    p0 / (1 + (r/rc)**2)
                ) 
            else:
                density.append(
                    ps / ((r/rs) * (1 + (r/rs))**2)
                )
        density = np.array(density)
    elif model == "nfw":
        rs = scale_length
        # Avoid division by zero
        radius_safe = np.where(radius == 0, 1e-10, radius)
        density = (
            p0 / ((radius_safe/rs) * (1 + (radius_safe/rs))**2)
        )
    
    return density

def mass_enclosed_function(scale_length, mass, radius, model, r1=None, rc=None, p0=None, ps=None):
    return 4 * np.pi * radius**2 * density_function(scale_length, mass, radius, model, r1, rc, p0, ps)

def counts_per_bin(scale_length, mass, radius, mass_per_particle, bin_width, model, r1=None, rc=None, p0=None, ps=None):
    return (
        mass_enclosed_function(scale_length, mass, radius, model, r1, rc, p0, ps)
        * (1 / mass_per_particle)
        * bin_width
    )

def kl_divergence(p, q):
    # Add small constant to avoid division by zero or log(0)
    epsilon = 1e-10
    p = np.array(p) + epsilon
    q = np.array(q) + epsilon
    
    # Normalize to ensure they sum to 1
    p = p / np.sum(p)
    q = q / np.sum(q)

    return np.sum(p * np.log(p / q))

def main():
    ####################################### Input parameters #######################################
    '''
    # Eric's Parameters
    radius_1 = 0.181216
    radius_2 = 0.182799
    mass_1 = 1.22251
    mass_2 = 0.0126171
    
    #Sid's Parameters
    radius_1 = 0.2
    radius_2 = 0.2
    mass_1 = 12.0
    mass_2 = 0.2
    '''
    
    radius_1 = 0.5
    radius_2 = 0.125
    mass_1 = 4.5
    mass_2 = 0.000136
    
    # Model types (keeping as plummer for now to maintain original functionality)
    model_l = "plummer"
    model_d = "cored"
    
    # Additional model parameters (None for plummer)
    r1_l = None
    rc_l = None
    r1_d = 0.7
    rc_d = 0.6

    # Number of particles
    nbodies_l = 0
    nbodies_d = 20000

    # Defining softening length array to test
    softening_length_array = np.logspace(-15, -1, 10) 
    
    # Number of timesteps and output frequence
    num_files = 3149
    output_frequency = 787

    # File name
    run_name = "victor"

    #################################################################################################

    my_path = os.path.abspath(os.path.dirname(__file__))

    # Defining scale radii and masses from ratio
    rscale_l = radius_1
    rscale_d = (radius_1 * (1 - radius_2)) / radius_2
    mass_l = mass_1
    mass_d = (mass_1 * (1 - mass_2)) / mass_2

    # Calculating timestep length in Gyr
    mass_enc_d = (
        mass_d * (rscale_l) ** 3 * ((rscale_l) ** 2 + (rscale_d) ** 2) ** (-3.0 / 2.0)
    )
    mass_enc_l = (
        mass_l * (rscale_d) ** 3 * ((rscale_l) ** 2 + (rscale_d) ** 2) ** (-3.0 / 2.0)
    )

    s1 = (rscale_l) ** 3 / (mass_enc_d + mass_l)
    s2 = (rscale_d) ** 3 / (mass_enc_l + mass_d)

    if s1 < s2:
        s = s1
    else:
        s = s2

    timestep = (1.0 / 100.0) * ((np.pi * (4.0 / 3.0) * s) ** (1.0 / 2.0))

    # Get extra variables if model is NFW or Cored 
    if model_l == "plummer":
        r200_l, p0_l, ps_l, bound_l = None, None, None, None
    else:
        r200_l, p0_l, ps_l, bound_l = set_vars(mass_l, rscale_l, model_l, r1_l, rc_l)
        mass_l = get_extra_nfw_mass(p0_l, ps_l, bound_l, model_l, rscale_l, r1_l, rc_l)
        print("mass_l: ", mass_l)
    if model_d == "plummer":
        r200_d, p0_d, ps_d, bound_d = None, None, None, None
    else:
        r200_d, p0_d, ps_d, bound_d = set_vars(mass_d, rscale_d, model_d, r1_d, rc_d)
        mass_d = get_extra_nfw_mass(p0_d, ps_d, bound_d, model_d, rscale_d, r1_d, rc_d)
        print("mass_d: ", mass_d)

    # Defining mass per particle
    if nbodies_l != 0:
        mass_per_particle_l = mass_l / nbodies_l
    else:
        mass_per_particle_l = 0
    if nbodies_d != 0:
        mass_per_particle_d = mass_d / nbodies_d
    else:
        mass_per_particle_d = 0

    
    # Defining bin width
    bin_width = min(rscale_l, rscale_d) / 5

    # Defining radius arrays (bin centers)
    if model_l == "plummer":
        search_radius_l = rscale_l * 4
    else:
        search_radius_l = rscale_l * 30
    if model_d == "plummer":
        search_radius_d = rscale_d * 4
    else:
        search_radius_d = rscale_d * 30
        
    radius_array_baryon = np.arange((rscale_l * 0.2), search_radius_l, bin_width)
    radius_array_dark = np.arange((rscale_d * 0.2), search_radius_d, bin_width)

    # Defining bin edges
    bin_edges_baryon = np.append(radius_array_baryon - bin_width/2, radius_array_baryon[-1] + bin_width/2)
    bin_edges_dark = np.append(radius_array_dark - bin_width/2, radius_array_dark[-1] + bin_width/2)

    # Defining array to store KL divergence values
    kl_divergence_array_baryon = []
    kl_divergence_array_dark = []

    # Defining array to store time values
    time_array = []

    # Defining array to store time values in Gyr
    for i in range(output_frequency - 1, num_files + 1, output_frequency): 
        time = i * timestep
        time_array.append(time)

    # Calculating theoretical density profiles
    if nbodies_l == 0:
        theoretical_baryon_density = np.zeros_like(radius_array_baryon)
    else:
        theoretical_baryon_density = counts_per_bin(
            rscale_l, mass_l, radius_array_baryon, mass_per_particle_l, bin_width,
            model_l, r1_l, rc_l, p0_l, ps_l
        )
    
    if nbodies_d == 0:
        theoretical_dark_density = np.zeros_like(radius_array_dark)
    else:
        theoretical_dark_density = counts_per_bin(
            rscale_d, mass_d, radius_array_dark, mass_per_particle_d, bin_width,
            model_d, r1_d, rc_d, p0_d, ps_d
        )

    # Calculating probability distributions for theoretical data
    if np.sum(theoretical_baryon_density) > 0:
        probability_distribution_baryon_theoretical = theoretical_baryon_density / np.sum(theoretical_baryon_density)
    else:
        probability_distribution_baryon_theoretical = np.zeros_like(theoretical_baryon_density)

    if np.sum(theoretical_dark_density) > 0:
        probability_distribution_dark_theoretical = theoretical_dark_density / np.sum(theoretical_dark_density)
    else:
        probability_distribution_dark_theoretical = np.zeros_like(theoretical_dark_density)

    for softening_length in softening_length_array:

        # Updating softening length in the Lua script
        with open(f"{my_path}/nbody/sample_workunits/for_developers.lua", "r") as file:
            lines = file.readlines()
            lines[218] = (f"      eps2        = {softening_length}, \n")

        with open(f"{my_path}/nbody/sample_workunits/for_developers.lua", "w") as file:
            file.writelines(lines)

        # Running the N-body simulation
        subprocess.run(["./run_nbody.sh"])

        # Defining array to store KL divergence values for this timestep
        kl_divergence_array_baryon_timestep = []
        kl_divergence_array_dark_timestep = []

        for i in range(output_frequency - 1, num_files + 1, output_frequency): 

            # Getting baryonic and dark matter data
            t = mwahpy.output_handler.read_output(f"{my_path}/build/bin/{i}") 

            light_x = []
            light_y = []
            light_z = []
            dark_x = []
            dark_y = []
            dark_z = []

            for j in range(len(t.x)):
                if t.typ[j] == 0:
                    light_x.append(t.x[j])
                    light_y.append(t.y[j])
                    light_z.append(t.z[j])
                else:
                    dark_x.append(t.x[j])
                    dark_y.append(t.y[j])
                    dark_z.append(t.z[j])

            light_x = np.array(light_x)
            light_y = np.array(light_y)
            light_z = np.array(light_z)
            dark_x = np.array(dark_x)
            dark_y = np.array(dark_y)
            dark_z = np.array(dark_z)

            light_r = np.sqrt(light_x**2 + light_y**2 + light_z**2)
            dark_r = np.sqrt(dark_x**2 + dark_y**2 + dark_z**2)

            # Calculate actual particle counts from simulation data using histograms
            if len(light_r) > 0:
                simulation_baryon_counts, _ = np.histogram(light_r, bins=bin_edges_baryon)
            else:
                simulation_baryon_counts = np.zeros(len(radius_array_baryon))
                
            if len(dark_r) > 0:
                simulation_dark_counts, _ = np.histogram(dark_r, bins=bin_edges_dark)
            else:
                simulation_dark_counts = np.zeros(len(radius_array_dark))

            # Calculating probability distributions for actual data
            if np.sum(simulation_baryon_counts) > 0:
                probability_distribution_baryon_simulation = simulation_baryon_counts / np.sum(simulation_baryon_counts)
            else:
                probability_distribution_baryon_simulation = np.zeros_like(simulation_baryon_counts)

            if np.sum(simulation_dark_counts) > 0:
                probability_distribution_dark_simulation = simulation_dark_counts / np.sum(simulation_dark_counts)
            else:
                probability_distribution_dark_simulation = np.zeros_like(simulation_dark_counts)

            # Calculating KL divergence between theoretical and actual distributions
            kl_baryon = kl_divergence(probability_distribution_baryon_simulation, probability_distribution_baryon_theoretical)
            kl_dark = kl_divergence(probability_distribution_dark_simulation, probability_distribution_dark_theoretical)

            # Storing KL divergence values for this timestep
            kl_divergence_array_baryon_timestep.append(kl_baryon)
            kl_divergence_array_dark_timestep.append(kl_dark)

        # Storing this timestep's KL divergence values
        kl_divergence_array_baryon.append(kl_divergence_array_baryon_timestep)
        kl_divergence_array_dark.append(kl_divergence_array_dark_timestep)

    # Convert to numpy arrays for easier manipulation
    kl_baryon_array = np.array(kl_divergence_array_baryon)
    kl_dark_array = np.array(kl_divergence_array_dark)
    softening_array = np.array(softening_length_array)
    softening_array = np.sqrt(softening_array)

    # Putting data into a pandas dataframe and saving to csv 
    df = pd.DataFrame({'softening_array': softening_array, 'kl_baryon': kl_baryon_array.tolist(), 'kl_dark': kl_dark_array.tolist()})
    df.to_csv(f"{my_path}/output/twenty_percent/{run_name}_soft_param_sweep_stability_test_timestep.csv", index=False)

    # Transpose the arrays so each row represents a timestep
    kl_baryon_array = np.array(kl_divergence_array_baryon).T
    kl_dark_array = np.array(kl_divergence_array_dark).T

    # Create a color map for the iterations
    colors = plt.cm.viridis(np.linspace(0, 1, len(kl_baryon_array)))

    # Plotting KL divergence values for Baryons
    plt.figure(figsize=(12, 6))  
    for i, color in enumerate(colors):
        plt.plot(softening_array, kl_baryon_array[i], 
                label=f'Baryons at {time_array[i]:.2f} Gyr', color=color)
    plt.axhline(y=0.0004944179238638797, color='teal', linestyle='--', label='Zeroth Timestep (Baryons)')
    plt.xscale('log')
    plt.xlabel('Softening Length')
    plt.ylabel('KL Divergence')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  
    plt.title('KL Divergence vs Softening Length for Baryons at Different Times')
    plt.tight_layout()  
    plt.savefig(f"{my_path}/output/twenty_percent/{run_name}_soft_param_sweep_stability_test_timestep_baryon.png", bbox_inches='tight') 
    plt.close()
    
    # Plotting KL divergence values for Dark Matter
    plt.figure(figsize=(12, 6))  
    for i, color in enumerate(colors):
        plt.plot(softening_array, kl_dark_array[i], 
                label=f'Dark Matter at {time_array[i]:.2f} Gyr', color=color)
    plt.axhline(y=0.0044291704720682374, color='pink', linestyle='--', label='Zeroth Timestep (Dark Matter)')
    plt.xscale('log')
    plt.xlabel('Softening Length')
    plt.ylabel('KL Divergence')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  
    plt.title('KL Divergence vs Softening Length for Dark Matter at Different Times')
    plt.tight_layout()  
    plt.savefig(f"{my_path}/output/twenty_percent/{run_name}_soft_param_sweep_stability_test_timestep_dark.png", bbox_inches='tight')  
    plt.close()

if __name__ == "__main__":
    main()