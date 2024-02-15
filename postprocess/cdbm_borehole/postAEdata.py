import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

def process_velocity_data(time, velocity):
    # Ensure uniform time intervals
    dt = time[1] - time[0]

    # Remove DC component
    velocity -= np.mean(velocity)

    # Perform FFT
    frequency = np.fft.fftfreq(len(time), dt)
    fft_result = np.fft.fft(velocity)
    amplitude_spectrum = np.abs(fft_result)

    return frequency, amplitude_spectrum, fft_result

def plot_frequency_spectrum(frequency, amplitude_spectrum, threshold=None):
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, amplitude_spectrum, label='Frequency Spectrum')
    
    if threshold:
        plt.axhline(y=threshold, color='r', linestyle='--', label='Threshold')

    plt.title('Frequency Spectrum Analysis')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.grid(True)
    plt.show()

def butter_highpass(cutoff, fs, order=4):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def highpass_filter(data, cutoff_frequency, fs, order=4):
    b, a = butter_highpass(cutoff_frequency, fs, order=order)
    filtered_data = filtfilt(b, a, data)
    return filtered_data

def plot_filtered_signal(time, original_velocity, filtered_signal, ptr_x, ptr_y, showfig=True, savefig=False):
    
    plt.figure(figsize=(10, 6))
    # plt.plot(time, original_velocity, label='Original Velocity', alpha=0.7)
    plt.plot(time, filtered_signal, linestyle='--', color='r')

    plt.title('Filtered Signal with Butterworth High-Pass Filter at location: '+str(ptr_x)+' , '+str(ptr_y))
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.grid(True)

    if savefig:
        plt.savefig('./outputs/img_filtered_signal'+str(i)+'.png')
    if showfig:
        plt.show()

    plt.close()

def plot_origin_signal(time, velocity, ptr_x, ptr_y, showfig=True, savefig=False):

    plt.figure(figsize=(10,6))
    plt.plot(time,velocity,'r-')
    plt.title('Original Signal at location: '+str(ptr_x)+' , '+str(ptr_y))
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.grid(True)

    if savefig:
        plt.savefig('./outputs/img_origin_signal'+str(i)+'.png')
    if showfig:
        plt.show()

    plt.close()

def plot_square_amplitude(time, filtered_velocity, ptr_x, ptr_y, showfig=True, savefig=False):

    plt.figure(figsize=(10,6))
    plt.stem(time,filtered_velocity*filtered_velocity,'r-')
    plt.title('Squared Amplitude of Filtered Velocity at location: '+str(ptr_x)+' , '+str(ptr_y))
    plt.xlabel('Time')
    plt.ylabel('Squared Amplitude Filtered Velocity')
    plt.grid(True)

    if savefig:
        plt.savefig('./outputs/img_squared_amplitude'+str(i)+'.png')
    if showfig:
        plt.show()
    
    plt.close()

# Example data
# time = np.linspace(0, 1, 1000)  # Replace with your actual time data
# velocity = np.sin(2 * np.pi * 10 * time) + 0.5 * np.random.randn(1000)  # Replace with your actual velocity data

# load points
given_coord_list = np.loadtxt("./outputs/given_coord_list.txt")

for i in range(np.shape(given_coord_list)[0]):

    #load files
    time = np.loadtxt("./outputs/list_timeseries"+str(i)+".txt")
    velocity = np.loadtxt("./outputs/list_velmag"+str(i)+".txt")

    #get point coordinates
    ptr_x = given_coord_list[i,0]
    ptr_y = given_coord_list[i,1]

    # Plot origin signal
    plot_origin_signal(time, velocity, ptr_x, ptr_y, showfig=True, savefig=True)

    # Process data
    frequency, amplitude_spectrum, _ = process_velocity_data(time, velocity)

    # Plot frequency spectrum
    plot_frequency_spectrum(frequency, amplitude_spectrum)

    # Set a cutoff frequency for the high-pass filter
    cutoff_frequency = 50000 #100000 # Adjust based on your data and requirements

    # Apply Butterworth high-pass filter
    filtered_velocity = highpass_filter(velocity, cutoff_frequency, fs=1/(time[1] - time[0]))

    # Plot the original and filtered signals
    plot_filtered_signal(time, velocity, filtered_velocity, ptr_x, ptr_y, showfig=True, savefig=True)

    # Plot square amplitude velocity
    plot_square_amplitude(time, filtered_velocity, ptr_x, ptr_y, showfig=True, savefig=True)