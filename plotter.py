import numpy as np
import mwahpy.output_handler
import matplotlib.pyplot as plt
import os.path
import imageio.v2 as imageio


def set_vars(mass, rscale, model, r1, rc):
    vol_pcrit = 0.568910904587397184785763397846734505212216314432372653620
    pcrit = 0.000679087369829744220469326744094105320596648627735869652
    r200 = (mass / vol_pcrit) ** (1.0 / 3.0)

    if model == "cored":
        D1 = r1 * (1 + r1 / rscale) ** 2 / (rscale + rscale * (r1 / rc) ** 2)
        D2 = (rscale) ** 3 * (
            np.log(1 + r200 / rscale)
            - np.log(1 + r1 / rscale)
            - r200 / (rscale + r200)
            + r1 / (rscale + r1)
        )
        D3 = (rc) ** 2 * (r1 - rc * np.arctan(r1 / rc))
        p0 = mass / (4 * np.pi * (D1 * D2 + D3))
        ps = p0 * D1
    elif model == "nfw":  # NFW
        c = r200 / rscale
        term = np.log(1.0 + c) - c / (1.0 + c)
        p0 = 200.0 * c**3 * pcrit / (3.0 * term)
        ps = None

    bound = 5.0 * r200

    return r200, p0, ps, bound


def get_extra_nfw_mass(p0, ps, bound, model, rscale, r1, rc):
    r = bound
    rs = rscale
    if model == "cored":
        C1 = 0
        C3 = C1 + 4 * np.pi * (
            ps * (rs**3) * (np.log(1 + r1 / rs) - r1 / (rs + r1))
            - p0 * (rc**2) * (r1 - rc * np.arctan(r1 / rc))
        )
        if r <= r1:
            m = 4.0 * np.pi * p0 * (rc**2) * (r - rc * np.arctan(r / rc)) - C1
        else:
            m = 4.0 * np.pi * (rs**3) * ps * (np.log(1 + r / rs) - r / (rs + r)) - C3
    elif model == "nfw":  # NFW
        m = 4.0 * np.pi * p0 * (rs**3) * (np.log((rs + r) / rs) - r / (rs + r))
    return m


def density_function(
    scale_length, mass, radius, model, r1=None, rc=None, p0=None, ps=None
):
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
                density.append(p0 / (1 + (r / rc) ** 2))
            else:
                density.append(ps / ((r / rs) * (1 + (r / rs)) ** 2))
        density = np.array(density)
    elif model == "nfw":
        rs = scale_length  # kpc
        # Avoid division by zero
        radius_safe = np.where(radius == 0, 1e-10, radius)
        density = p0 / ((radius_safe / rs) * (1 + (radius_safe / rs)) ** 2)
    elif model == "hernquist":
        density = (
            1
            / (2 * np.pi)
            * mass
            * scale_length
            / (radius * (radius + scale_length) ** 3)
        )

    return density


def mass_enclosed_function(
    scale_length, mass, radius, model, r1=None, rc=None, p0=None, ps=None
):
    return (
        4
        * np.pi
        * radius**2
        * density_function(scale_length, mass, radius, model, r1, rc, p0, ps)
    )


def counts_per_histogram_bin(
    scale_length,
    mass,
    radius,
    mass_per_particle,
    bin_width,
    model,
    r1=None,
    rc=None,
    p0=None,
    ps=None,
):
    return (
        mass_enclosed_function(scale_length, mass, radius, model, r1, rc, p0, ps)
        * (1 / mass_per_particle)
        * bin_width
    )


def single_scatter_plot(
    x,
    y,
    z,
    time,
    my_path,
    save_dir,
    type,
    index,
    parameter_name,
    run_type,
    softening_parameter,
    model_l,
    model_d,
    rscale_l,
    rscale_d,
):
    if type == "baryonic":
        color = "blue"
        model = model_l
        rscale = rscale_l
    elif type == "dark":
        color = "red"
        model = model_d
        rscale = rscale_d
    else:
        color = "purple"
        if model_l == "plummer" and model_d == "plummer":
            model = model_l
        else:
            model = model_l if model_l != "plummer" else model_d
        rscale = max(rscale_l, rscale_d)

    if model == "plummer":
        axis_lim = 10 * rscale
    else:
        axis_lim = 100 * rscale

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(
        f"{type.capitalize()} Matter: {parameter_name}'s Parameters with {run_type.replace('_', ' ').title()} and {softening_parameter.replace('_', ' ').title()}"
    )
    fig.add_subplot(2, 2, 1)
    plt.scatter(x, y, marker=".", s=0.1, alpha=0.5, c=color)
    plt.title("X vs Y")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("x (kpc)")
    plt.ylabel("y (kpc)")
    fig.add_subplot(2, 2, 2)
    plt.scatter(y, z, marker=".", s=0.1, alpha=0.5, c=color)
    plt.title("Y vs Z")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("y (kpc)")
    plt.ylabel("z (kpc)")
    fig.add_subplot(2, 2, 3)
    plt.scatter(x, z, marker=".", s=0.1, alpha=0.5, c=color)
    plt.title("X vs Z")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("x (kpc)")
    plt.ylabel("z (kpc)")
    fig.add_subplot(2, 2, 4)
    plt.text(0.5, 0.5, f"Time: {time} Gyr", fontsize=12, ha="center")
    plt.axis("off")
    plt.subplots_adjust(hspace=0.5)
    fig_save_name = f"{parameter_name}_{run_type}_{softening_parameter}_{type}_matter_cartisean_position"
    fig.savefig(f"{my_path}/{save_dir}/scatter/{fig_save_name}_{str(index)}_s.png")
    plt.close()


def seperated_scatter_plot(
    light_x,
    light_y,
    light_z,
    dark_x,
    dark_y,
    dark_z,
    time,
    my_path,
    save_dir,
    index,
    parameter_name,
    run_type,
    softening_parameter,
    model_l,
    model_d,
    rscale_l,
    rscale_d,
):
    if model_l == "plummer" and model_d == "plummer":
        axis_lim = 10 * max(rscale_l, rscale_d)
    else:
        axis_lim = 100 * max(rscale_l, rscale_d)

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(
        f"Seperated Matter: {parameter_name}'s Parameters with {run_type.replace('_', ' ').title()} and {softening_parameter.replace('_', ' ').title()}"
    )
    fig.add_subplot(2, 2, 1)
    plt.scatter(dark_x, dark_y, marker=".", s=0.1, alpha=0.5, c="red")
    plt.scatter(light_x, light_y, marker=".", s=0.1, alpha=0.5, c="blue")
    plt.title("X vs Y")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("x (kpc)")
    plt.ylabel("y (kpc)")
    fig.add_subplot(2, 2, 2)
    plt.scatter(dark_y, dark_z, marker=".", s=0.5, alpha=0.1, c="red")
    plt.scatter(light_y, light_z, marker=".", s=0.5, alpha=0.1, c="blue")
    plt.title("Y vs Z")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("y (kpc)")
    plt.ylabel("z (kpc)")
    fig.add_subplot(2, 2, 3)
    plt.scatter(dark_x, dark_z, marker=".", s=0.5, alpha=0.1, c="red")
    plt.scatter(light_x, light_z, marker=".", s=0.5, alpha=0.1, c="blue")
    plt.title("X vs Z")
    plt.axis([-axis_lim, axis_lim, -axis_lim, axis_lim])
    plt.xlabel("x (kpc)")
    plt.ylabel("z (kpc)")
    fig.add_subplot(2, 2, 4)
    plt.text(0.5, 0.5, f"Time: {time} Gyr", fontsize=12, ha="center")
    plt.axis("off")
    plt.subplots_adjust(hspace=0.5)
    fig_save_name = f"{parameter_name}_{run_type}_{softening_parameter}_seperated_matter_cartisean_position"
    fig.savefig(f"{my_path}/{save_dir}/scatter/{fig_save_name}_{str(index)}_s.png")
    plt.close()


def density_profile_hist(
    rscale_l,
    rscale_d,
    mass_l,
    mass_d,
    light_r,
    dark_r,
    model_l,
    model_d,
    combined_r,
    radius_array_hist,
    radius_array_theoretical,
    mass_per_particle_l,
    mass_per_particle_d,
    nbodies_l,
    nbodies_d,
    time,
    my_path,
    save_dir,
    parameter_name,
    run_type,
    softening_parameter,
    index,
    bin_width,
    r1_l=None,
    rc_l=None,
    p0_l=None,
    ps_l=None,
    r1_d=None,
    rc_d=None,
    p0_d=None,
    ps_d=None,
):
    # Calculating theoretical density profiles
    if nbodies_l != 0:
        theoretical_baryon_density = counts_per_histogram_bin(
            rscale_l,
            mass_l,
            radius_array_theoretical,
            mass_per_particle_l,
            bin_width,
            model_l,
            r1_l,
            rc_l,
            p0_l,
            ps_l,
        )
    else:
        theoretical_baryon_density = np.zeros(len(radius_array_theoretical))
    if nbodies_d != 0:
        theoretical_dark_density = counts_per_histogram_bin(
            rscale_d,
            mass_d,
            radius_array_theoretical,
            mass_per_particle_d,
            bin_width,
            model_d,
            r1_d,
            rc_d,
            p0_d,
            ps_d,
        )
    else:
        theoretical_dark_density = np.zeros(len(radius_array_theoretical))
    theoretical_combined_density = theoretical_dark_density + theoretical_baryon_density

    # Calculating axis limits
    if model_l == "plummer" or model_l == "hernquist":
        baryon_xlim = 10 * rscale_l
    else:
        baryon_xlim = 100 * rscale_l

    if model_d == "plummer" or model_d == "hernquist":
        dark_xlim = 10 * rscale_d
    else:
        dark_xlim = 100 * rscale_d
    combined_xlim = max(baryon_xlim, dark_xlim)
    baryon_ylim = 1.1 * (max(theoretical_baryon_density))
    dark_ylim = 1.1 * (max(theoretical_dark_density))
    combined_ylim = 1.1 * (max(theoretical_combined_density))

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(
        f"Density Profile: {parameter_name}'s Parameters with {run_type.replace('_', ' ').title()} and {softening_parameter.replace('_', ' ').title()}"
    )
    fig.add_subplot(2, 2, 1)
    plt.plot([], [], " ", label=f"{nbodies_l} Particles")
    plt.plot(
        radius_array_theoretical,
        theoretical_baryon_density,
        label="Theoretical (Baryonic)",
        color="blue",
    )
    plt.hist(
        light_r,
        bins=np.arange(
            min(radius_array_hist), max(radius_array_hist) + bin_width / 2, bin_width
        ),
        label="Simulation (Baryonic)",
        color="teal",
    )
    plt.legend()
    plt.xlim(0, baryon_xlim)
    plt.ylim(0, baryon_ylim)
    plt.ylabel("N")
    plt.title("Baryonic Matter Density Profile")
    fig.add_subplot(2, 2, 2)
    plt.plot([], [], " ", label=f"{nbodies_d} Particles")
    plt.plot(
        radius_array_theoretical,
        theoretical_dark_density,
        label="Theoretical (Dark)",
        color="red",
    )
    plt.hist(
        dark_r,
        bins=np.arange(
            min(radius_array_hist), max(radius_array_hist) + bin_width / 2, bin_width
        ),
        label="Simulation (Dark)",
        color="pink",
    )
    plt.xlim(0, dark_xlim)
    plt.ylim(0, dark_ylim)
    plt.legend()
    plt.title("Dark Matter Density Profile")
    fig.add_subplot(2, 2, 3)
    plt.plot([], [], " ", label=f"{nbodies_l + nbodies_d} Particles")
    plt.plot(
        radius_array_theoretical,
        theoretical_combined_density,
        label="Theoretical (Combined)",
        color="purple",
    )
    plt.hist(
        combined_r,
        bins=np.arange(
            min(radius_array_hist), max(radius_array_hist) + bin_width / 2, bin_width
        ),
        label="Simulation (Combined)",
        color="#D8BFD8",
    )
    plt.xlim(0, combined_xlim)
    plt.ylim(0, combined_ylim)
    plt.ylabel("N")
    plt.legend()
    plt.xlabel("r (kpc)")
    plt.title("Combined Density Profile")
    fig.add_subplot(2, 2, 4)
    plt.plot(
        radius_array_theoretical,
        theoretical_combined_density,
        label="Theoretical (Combined)",
        color="purple",
    )
    plt.plot(
        radius_array_theoretical,
        theoretical_baryon_density,
        label="Theoretical (Baryonic)",
        color="blue",
    )
    plt.plot(
        radius_array_theoretical,
        theoretical_dark_density,
        label="Theoretical (Dark)",
        color="red",
    )
    if max(theoretical_baryon_density) > max(theoretical_dark_density):
        plt.hist(
            combined_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Combined)",
            color="#D8BFD8",
        )
        plt.hist(
            light_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Baryonic)",
            color="teal",
        )
        plt.hist(
            dark_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Dark)",
            color="pink",
        )
    else:
        plt.hist(
            combined_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Combined)",
            color="#D8BFD8",
        )
        plt.hist(
            dark_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Dark)",
            color="pink",
        )
        plt.hist(
            light_r,
            bins=np.arange(
                min(radius_array_hist),
                max(radius_array_hist) + bin_width / 2,
                bin_width,
            ),
            label="Simulation (Baryonic)",
            color="teal",
        )
    plt.legend()
    plt.xlim(0, max([baryon_xlim, dark_xlim, combined_xlim]))
    plt.ylim(0, max([baryon_ylim, dark_ylim, combined_ylim]))
    plt.xlabel("r (kpc)")
    plt.title("Overlayed Density Profiles")
    plt.subplots_adjust(hspace=0.5)
    fig.text(0.5, 0.5, f"Time: {time} Gyr", fontsize=12, ha="center", va="center")
    fig_save_name = f"{parameter_name}_{run_type}_{softening_parameter}_density_profile"
    fig.savefig(f"{my_path}/{save_dir}/hist/{fig_save_name}_{str(index)}_h.png")
    plt.close()


def make_movie(
    my_path,
    save_dir,
    parameter_name,
    start_index,
    end_index,
    output_frequency,
    graph_type,
    matter_type,
    run_type,
    softening_parameter,
):
    if graph_type == "scatter":
        filenames = [
            f"{my_path}/{save_dir}/scatter/{parameter_name}_{run_type}_{softening_parameter}_{matter_type}_matter_cartisean_position_{str(index)}_s.png"
            for index in range(start_index, end_index, output_frequency)
        ]
        imageio.mimsave(
            f"{save_dir}/gif/{parameter_name}_{run_type}_{softening_parameter}_{matter_type}_matter_animated_scatter_plot.gif",
            [imageio.imread(filename) for filename in filenames],
            duration=1,
        )
    else:
        filenames = [
            f"{my_path}/{save_dir}/hist/{parameter_name}_{run_type}_{softening_parameter}_density_profile_{str(index)}_h.png"
            for index in range(start_index, end_index, output_frequency)
        ]
        imageio.mimsave(
            f"{save_dir}/gif/{parameter_name}_{run_type}_{softening_parameter}_animated_density_plot.gif",
            [imageio.imread(filename) for filename in filenames],
            duration=1,
        )


def main():
    ####################################### Input parameters #######################################
    num_files = 7699
    output_frequency = 100

    """
    # Eric's Parameters
    radius_1 = 0.181216
    radius_2 = 0.182799
    mass_1 = 1.22251
    mass_2 = 0.0126171
    """
    # Sidd's Parameters
    radius_1 = 0.2
    radius_2 = 0.2
    mass_1 = 12.0
    mass_2 = 0.2
    """
    # Fake LeoT Parameters
    radius_1 = 0.1
    radius_2 = 0.3
    mass_1 = 0.45
    mass_2 = 0.1
    
    # SGR Parameters
    radius_1 = 1.0
    radius_2 = 0.1
    mass_1 = 450.0
    mass_2 = 0.01
    
    # Victor's Parameters
    radius_1 = 0.2
    radius_2 = 0.054054
    mass_1 = 4.5
    mass_2 = 0.000136345
    """
    # Cored Parameters (not used for any other model)
    r1_l = 0.7
    rc_l = 0.6
    r1_d = 0.7
    rc_d = 0.6

    nbodies_l = 100
    nbodies_d = 100
    model_l = "plummer"
    model_d = "plummer"

    parameter_name = "Test"
    run_type = "small"
    extra_parameter = ""
    if extra_parameter == "":
        save_dir = f"{parameter_name}/{run_type}"
    else:
        save_dir = f"{parameter_name}/{run_type}/{extra_parameter}"
    #################################################################################################

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

    my_path = os.path.abspath(os.path.dirname(__file__))

    # Set bin width
    bin_width = min(rscale_l, rscale_d) / 5

    # Set radius array for theoretical data and simulation data histogram
    if model_l == "plummer" and model_d == "plummer":
        radius_array_theoretical = np.arange(1e-6, 10 * max(rscale_l, rscale_d), 0.001)
        radius_array_hist = np.concatenate(
            (
                np.array([0]),
                np.arange(bin_width / 2, 10 * max(rscale_l, rscale_d), bin_width),
            )
        )
    else:
        radius_array_theoretical = np.arange(1e-6, 1000 * rscale_l, 0.001)
        radius_array_hist = np.concatenate(
            (
                np.array([0]),
                np.arange(bin_width / 2, 150 * max(rscale_l, rscale_d), bin_width),
            )
        )

    # Get extra variables if model is NFW or Cored
    if model_l == "plummer" or model_l == "hernquist":
        r200_l, p0_l, ps_l, bound_l = None, None, None, None
    else:
        r200_l, p0_l, ps_l, bound_l = set_vars(mass_l, rscale_l, model_l, r1_l, rc_l)
        mass_l = get_extra_nfw_mass(p0_l, ps_l, bound_l, model_l, rscale_l, r1_l, rc_l)
        print("mass_l: ", mass_l)
    if model_d == "plummer" or model_d == "hernquist":
        r200_d, p0_d, ps_d, bound_d = None, None, None, None
    else:
        r200_d, p0_d, ps_d, bound_d = set_vars(mass_d, rscale_d, model_d, r1_d, rc_d)
        mass_d = get_extra_nfw_mass(p0_d, ps_d, bound_d, model_d, rscale_d, r1_d, rc_d)
        print("mass_d: ", mass_d)

    # Calculating mass per particle
    if nbodies_l != 0:
        mass_per_particle_l = mass_l / nbodies_l
    else:
        mass_per_particle_l = 0
    if nbodies_d != 0:
        mass_per_particle_d = mass_d / nbodies_d
    else:
        mass_per_particle_d = 0

    # Check if initial file exists and determine starting index
    initial_file_path = f"{save_dir}/data/initial.out"
    if os.path.exists(initial_file_path):
        start_index = -1
        initial_file_processed = True
    else:
        start_index = output_frequency - 1
        initial_file_processed = False

    for i in range(start_index, num_files + 1, output_frequency):
        # Handle special case for initial file
        if i == -1:
            filename = "initial.out"
        else:
            filename = str(i)

        t = mwahpy.output_handler.read_output(f"{save_dir}/data/{filename}")
        # For initial file, time is 0, otherwise use the timestep calculation
        time = 0 if i == -1 else i * timestep
        
        # Convert to numpy arrays and separate by type using boolean indexing
        x = np.array(t.x)
        y = np.array(t.y)
        z = np.array(t.z)
        typ = np.array(t.typ)
        
        # Create boolean masks for light (type 0) and dark (type 1) matter
        light_mask = typ == 0
        dark_mask = typ == 1
        
        # Extract light and dark matter data using boolean indexing
        light_x = x[light_mask]
        light_y = y[light_mask]
        light_z = z[light_mask]
        dark_x = x[dark_mask]
        dark_y = y[dark_mask]
        dark_z = z[dark_mask]
        
        # Combined data (all particles)
        combined_x = x
        combined_y = y
        combined_z = z

        # Making scatter and histogram plots
        if nbodies_l != 0:
            single_scatter_plot(
                light_x,
                light_y,
                light_z,
                time,
                my_path,
                save_dir,
                "baryonic",
                i,
                parameter_name,
                run_type,
                extra_parameter,
                model_l,
                model_d,
                rscale_l,
                rscale_d,
            )
        if nbodies_d != 0:
            single_scatter_plot(
                dark_x,
                dark_y,
                dark_z,
                time,
                my_path,
                save_dir,
                "dark",
                i,
                parameter_name,
                run_type,
                extra_parameter,
                model_l,
                model_d,
                rscale_l,
                rscale_d,
            )
        if nbodies_l != 0 and nbodies_d != 0:
            single_scatter_plot(
                combined_x,
                combined_y,
                combined_z,
                time,
                my_path,
                save_dir,
                "combined",
                i,
                parameter_name,
                run_type,
                extra_parameter,
                model_l,
                model_d,
                rscale_l,
                rscale_d,
            )
            seperated_scatter_plot(
                light_x,
                light_y,
                light_z,
                dark_x,
                dark_y,
                dark_z,
                time,
                my_path,
                save_dir,
                i,
                parameter_name,
                run_type,
                extra_parameter,
                model_l,
                model_d,
                rscale_l,
                rscale_d,
            )

        # Calculate radial distances (arrays are already numpy arrays)
        light_r = np.sqrt(light_x**2 + light_y**2 + light_z**2)
        dark_r = np.sqrt(dark_x**2 + dark_y**2 + dark_z**2)
        combined_r = np.concatenate((light_r, dark_r))

        density_profile_hist(
            rscale_l,
            rscale_d,
            mass_l,
            mass_d,
            light_r,
            dark_r,
            model_l,
            model_d,
            combined_r,
            radius_array_hist,
            radius_array_theoretical,
            mass_per_particle_l,
            mass_per_particle_d,
            nbodies_l,
            nbodies_d,
            time,
            my_path,
            save_dir,
            parameter_name,
            run_type,
            extra_parameter,
            i,
            bin_width,
            r1_l,
            rc_l,
            p0_l,
            ps_l,
            r1_d,
            rc_d,
            p0_d,
            ps_d,
        )

    # Making gif
    if nbodies_l == 0:
        matter_type_list = ["dark"]
    elif nbodies_d == 0:
        matter_type_list = ["baryonic"]
    else:
        matter_type_list = ["baryonic", "dark", "combined", "seperated"]

    # Adjust start index based on whether initial file was processed
    movie_start_index = -1 if initial_file_processed else start_index

    for matter_type in matter_type_list:
        make_movie(
            my_path,
            save_dir,
            parameter_name,
            movie_start_index,  # Use adjusted start index
            num_files + 1,
            output_frequency,
            "scatter",
            matter_type,
            run_type,
            extra_parameter,
        )

    make_movie(
        my_path,
        save_dir,
        parameter_name,
        movie_start_index,  # Use adjusted start index
        num_files + 1,
        output_frequency,
        "hist",
        "None",
        run_type,
        extra_parameter,
    )


if __name__ == "__main__":
    main()
