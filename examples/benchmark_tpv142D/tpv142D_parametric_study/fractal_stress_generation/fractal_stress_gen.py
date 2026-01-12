import numpy as np

def generate_fractal(length, D, target_mean, target_std, seed):
    """
    Generate a fractal shear stress profile using an inverse FFT.
    
    Parameters:
        length (int): Number of points in the shear stress profile.
        D (float): Fractal dimension of the profile.
        target_mean (float): Desired mean value of the output profile.
        target_std (float): Desired standard deviation of the output profile.
        seed (int): Seed for random number generation.
    
    Returns:
        np.ndarray: A 1D array representing the shear stress profile.
    """
    # Compute the Hurst exponent
    H = 2 - D

    # Set the random seed for reproducibility
    np.random.seed(seed)
    
    # Prepare the spectrum array.
    # For a real FFT inverse (irfft), the length of the complex spectrum is length//2 + 1.
    spectrum = np.zeros(length // 2 + 1, dtype=np.complex128)
    
    # Frequency domain: calculate power spectrum and assign a random phase to each component.
    for i in range(0, length // 2 + 1):
        # Avoid division by zero for the zero frequency term.
        freq = 1.0 if i == 0 else i / length
        power = 1.0 / (freq ** (H + 0.5))
        
        # Generate a random phase uniformly in [0, 2*pi)
        phase = np.random.uniform(0.0, 2 * np.pi)
        
        # Create the complex number for the spectrum with given power and phase.
        spectrum[i] = power * (np.cos(phase) + 1j * np.sin(phase))
    
    # Perform the inverse FFT to obtain the shear stress profile.
    shear_stress = np.fft.irfft(spectrum, n=length)
    
    # Normalize by the length (mimicking FFTW's scaling behavior)
    shear_stress /= length
    
    # Normalize the shear stress profile to have the target mean and standard deviation.
    current_mean = np.mean(shear_stress)
    current_std = np.std(shear_stress)
    shear_stress = (shear_stress - current_mean) / current_std * target_std + target_mean
    
    return shear_stress

# Example usage:
if __name__ == "__main__":

    #main fault point number: 560
    #branch fault point number: 238

    fault_name = "main_fault"  # or "branch_fault"
    number_of_files = 100  # Number of files to process

    D = 1.3              # Fractal dimension for the profile
    target_mean = 70e6   # Desired mean
    target_std = 5e6     # Desired standard deviation

    for i in range(number_of_files):

        seed = np.random.randint(0, 2**32)  # Random seed for reproducibility

        if fault_name == "main_fault":

            coords_file = "../data_used_extract_qp_locs/main_fault_qp_coords.csv"
            output_filename = "main_fault_fractal_stress_"+str(i)+".csv"

            numofptr = 560        # Number of points in the profile

            # Generate fractal shear stress profile
            profile = generate_fractal(numofptr, D, target_mean, target_std, seed)

            # read the coordinates of the main fault
            coords = np.loadtxt(coords_file, delimiter=',', skiprows=1)  # Assuming the file has a header

            # stack the coordinates and profile
            data = np.column_stack((coords, profile))
            # save the data to a CSV file
            np.savetxt(output_filename, data, delimiter=',', header='x,y,z,shear_stress', comments='')
            print(f"Fractal shear stress profile saved to {output_filename}")
        
        elif fault_name == "branch_fault":

            coords_file = "../data_used_extract_qp_locs/branch_fault_qp_coords.csv"
            output_filename = "branch_fault_fractal_stress_"+str(i)+".csv"
            
            numofptr = 238        # Number of points in the profile

            # Generate fractal shear stress profile
            profile = generate_fractal(numofptr, D, target_mean, target_std, seed)

            # read the coordinates of the main fault
            coords = np.loadtxt(coords_file, delimiter=',', skiprows=1)  # Assuming the file has a header

            # stack the coordinates and profile
            data = np.column_stack((coords, profile))
            # save the data to a CSV file
            np.savetxt(output_filename, data, delimiter=',', header='x,y,z,shear_stress', comments='')
            print(f"Fractal shear stress profile saved to {output_filename}")
        
        # # Optionally, visualize the profile (requires matplotlib)
        # import matplotlib.pyplot as plt
        # plt.plot(profile)
        # plt.title("Fractal Shear Stress Profile")
        # plt.xlabel("Index")
        # plt.ylabel("Shear Stress")
        # plt.show()