from paul import *
from plotting import plot_parameters, plot_simulation_results
from checks import run_checks


def run_generation(generation):
    get_couples(generation)
    get_generation_embryos(generation)
    select_generation_embryos(generation)


def save_results():
    save_generation_scores()
    save_intergenerational_scores()
    save_freq_vs_weight()
    store_ukb_pca_map()
    store_ukb_pca_array()
    store_homozygosity()  # 136 minutes


def main():
    
    # Store the preparatory data
    store_initial_generation()

    # Plot the parameters governing the simulation
    plot_parameters()

    # Run the simulation over generations
    for generation in range(N_GENERATIONS):
        run_generation(generation)

    # Save the results
    save_results()

    # Plot the results
    plot_simulation_results()

    run_checks()



if __name__ == "__main__":
    main()



