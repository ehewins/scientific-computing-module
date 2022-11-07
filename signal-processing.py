import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from scipy.signal import square, sawtooth

# NOTE TO EXAMINER: Every feature of this program seems to work except for the windowing. As soon as you change the time or frequency windowing, you start to get buggy behaviour.
# The inverse fourier transform is dependent on the windowing, and only produces the correct graph when the windoing is set to include the full range of data.
# When you change the slider for the total time, you should also move the slider for the time and frequency window upper limits to their maximum values, otherwise the new data doesn't get included in the calculations.
# I went to Kasper's office to ask for help on Monday at 12:00, but we couln't fix the probelms, but he said the other aspects were solid. I would like to submit this code again after the deadline with fixes made, but I've run out of time for now.

def generate_signal(type, A, f, t, phase):
    '''Function to generate a signal of type Sine, Square or Sawtooth when given the amplitude, frequency [Hz], an array of time values [s], and the phase [rad]'''
    if type == 'Sine':
        return A*np.sin(2*np.pi*f*t + phase)
    elif type == 'Square':
        return A*square(2*np.pi*f*t + phase)
    elif type == 'Sawtooth':
        return A*sawtooth(2*np.pi*f*t + phase)
    else:
        print("Supplied signal type was invalid")
        return 0

def update_plots(val):
    '''Funcion which re-performs all the signal processing calculations and re-plots the graphs when any of the parameters are changed by the user'''
    # Obtain values of all the varibles from the positions of the sliders
    N = samples_handle.val
    T = ttot_handle.val
    sig_f1 = frequency1_handle.val
    A1 = amp1_handle.val
    phase1 = phase1_handle.val
    sig_f2 = frequency2_handle.val
    A2 = amp2_handle.val
    phase2 = phase2_handle.val
    signal_type = waveform_handle.value_selected
    noise_amp = noise_handle.val
    t0 = t0_handle.val
    t1 = t1_handle.val
    f0 = f0_handle.val
    f1 = f1_handle.val
    dt = T/N

    # generate an array of the time values
    times = np.arange(0, T, dt)
    # generate the signal using supplied values
    signal = generate_signal(signal_type, A1, sig_f1, times, phase1) + generate_signal(signal_type, A2, sig_f2, times, phase2)
    if noise_amp != 0: # only calculate noise if the noise slider is set to a non-zero value
        signal += 2*noise_amp*np.random.random_sample(N)-noise_amp
    signal_handle.set_data(times, signal) # Plot the data
    signal_ax.set_xlim(0,T)

    # Make sure the time windowing sliders & their values can't pass through each other.
    t0_handle.valmax = t1
    t1_handle.valmin = t0
    t1_handle.valmax = T
    t0_handle.ax.set_xlim(0,T) # Also reset the slider limits to account for changing T
    t1_handle.ax.set_xlim(0,T)
    # Give a visual representation of the windowing limits on the graph
    t0_line.set_xdata([t0,t0])
    t1_line.set_xdata([t1,t1])
    # Limit the data included in the future calculations based on the windowing
    signal = signal[(times >= t0) & (times <= t1)]
    times = times[(times >= t0) & (times <= t1)]
    # Redefine new N and T before calculating the frequencies
    dT = T/N # spacicing is the same as prior to windowing, so should be calculated with the old values before they're redefined on the next two lines
    N = np.size(times)
    T = N*dT

    # produce an array of the possible frequencies in the signal
    frequencies = (1/T * np.arange(N))
    for k in range(N): frequencies[k] -= N if k > N/2 else 0 # converts the second half of the array to their negative frequency counterparts
    nyquist = 0.5 * N/T

    # Perform the Fourier transform and plot the power spectrum
    FT_sig = np.fft.fft(signal, N)
    power_spectrum = (np.real(FT_sig)**2 + np.imag(FT_sig)**2)/(N**2)
    power_handle.set_data(frequencies[frequencies>=0], power_spectrum[frequencies>=0])
    power_ax.set_xlim(0,nyquist)
    power_ax.set_ylim(0,np.max(power_spectrum)*1.1)

    # Make sure the frequency windowing sliders & their values can't pass through each other.
    f0_handle.valmax = f1
    f1_handle.valmin = f0
    f1_handle.valmax = nyquist
    f0_handle.ax.set_xlim(0,nyquist) # Also reset the slider limits
    f1_handle.ax.set_xlim(0,nyquist)
    # Give a visual representation of the windowing limits on the graph
    f0_line.set_xdata([f0,f0])
    f1_line.set_xdata([f1,f1])
    # Limit the data included in the future calculations based on the windowing
    # should only keep the positive and negative components whose magnitudes are within the windowing limits
    FT_sig = FT_sig[(abs(frequencies) >= f0) & (abs(frequencies) <= f1)]
    frequencies = frequencies[(abs(frequencies) >= f0) & (abs(frequencies) <= f1)]

    # Inverse fourier transform from windowed frequency data. Doesn't work if the windowing has been changed from default values unfortunately
    IFT_sig = np.fft.ifft(FT_sig, N)
    inverse_handle.set_data(times, np.real(IFT_sig))
    inverse_ax.set_xlim(t0,t1)

    # redraw the plots
    fig.canvas.draw_idle()

def close_callback(event):
    '''Function used to close the figure window and exit the program, called via a button widget'''
    plt.close("all")
    exit()

# Initial values fore core variables:
T = 1 # The total time frame of the signal in seconds
sig_f = 10 # The frequency of the signal in Hertz
N = 100 # The number of sampling points
A = 1 # The amplitude of the signal in arbitrary units
A_max = 3 # The maximum amplitude of the signal, as well as of the noise
dt = T/N # spacing between time values

times = np.arange(0, T, dt) # Array of time values
signal = generate_signal('Sine', A, sig_f, times, 0) # The signal

# Generating the frequencies corresoponding to the Fourier transformed signal
frequencies = (1/T * np.arange(N))
for k in range(N): frequencies[k] -= N if k > N/2 else 0
nyquist = 0.5 * N/T

FT_sig = np.fft.fft(signal, N) # Fourier transform of sigal
power_spectrum = (np.real(FT_sig)**2 + np.imag(FT_sig)**2)/(N**2) # Generate the power spectrum
# Producing the inverse Fourier transform from the Fourier transform. Should be idential to the original signal until you start changing the windowing.
IFT_sig = np.fft.ifft(FT_sig)

fig = plt.figure(1,figsize=(16,9)) # Setting up the figure window

# Setup for the signal graph
signal_ax = plt.axes([0.04, 0.5, 0.28, 0.4])
signal_ax.set_title("Signal")
signal_ax.set_xlabel("Time, s")
signal_ax.set_ylabel("Amplitude")
signal_ax.set_xlim(0,T)
signal_ax.set_ylim(-A_max-0.1, A_max+0.1)
signal_handle, = signal_ax.plot(times, signal, 'r.-')

# Setup for the power spectrum graph
power_ax = plt.axes([0.365, 0.5, 0.28, 0.4])
power_ax.set_title("Power spectrum")
power_ax.set_xlabel("Frequency, Hz")
power_ax.set_ylabel("Amplitude$^2$")
power_ax.set_xlim(0,nyquist)
# power_ax.set_xlim(min(frequencies),max(frequencies))
power_ax.set_ylim(0,np.max(power_spectrum)*1.1)
power_handle, = power_ax.plot(frequencies[frequencies>=0], power_spectrum[frequencies>=0], 'b.-')

# Setup for the inverse fourier transform graph
inverse_ax = plt.axes([0.69, 0.5, 0.28, 0.4])
inverse_ax.set_title("Inverse FT")
inverse_ax.set_xlabel("Time, s")
inverse_ax.set_ylabel("Amplitude")
inverse_ax.set_xlim(0,T)
inverse_ax.set_ylim(-A_max-0.1, A_max+0.1)
inverse_handle, = inverse_ax.plot(times, np.real(IFT_sig), 'r.-')

# Setup for the close button
close_ax = plt.axes([0.79,0.35,0.08,0.06])
close_handle = widgets.Button(close_ax, 'Exit')
close_handle.on_clicked(close_callback)

# Setup for the num. samples slider
samples_ax = plt.axes([0.08, 0.42, 0.22, 0.02])
samples_handle = widgets.Slider(samples_ax, 'Num. samples', 10, 300, valinit=100, valstep=1)
samples_handle.on_changed(update_plots)

# Setup for the total time slider
ttot_ax = plt.axes([0.08, 0.39, 0.22, 0.02])
ttot_handle = widgets.Slider(ttot_ax, 'Total time [s]', 1, 4, valinit=1)
ttot_handle.on_changed(update_plots)

# Setup for the frequency slider of the first signal component
frequency1_ax = plt.axes([0.08, 0.34, 0.22, 0.02])
frequency1_handle = widgets.Slider(frequency1_ax, 'Frequency 1 [Hz]', 0, 20, valinit=10)
frequency1_handle.on_changed(update_plots)

# Setup for the amplitude slider of the first signal component
amp1_ax = plt.axes([0.08, 0.31, 0.22, 0.02])
amp1_handle = widgets.Slider(amp1_ax, 'Amplitude 1', 0, A_max, valinit=1)
amp1_handle.on_changed(update_plots)

# Setup for the phase slider of the first signal component
phase1_ax = plt.axes([0.08, 0.28, 0.22, 0.02])
phase1_handle = widgets.Slider(phase1_ax, 'Phase 1 [rad]', 0., 2*np.pi, valinit=0)
phase1_handle.on_changed(update_plots)

# Setup for the frequency slider of the second signal component
frequency2_ax = plt.axes([0.08, 0.23, 0.22, 0.02])
frequency2_handle = widgets.Slider(frequency2_ax, 'Frequency 2 [Hz]', 0, 20, valinit=0)
frequency2_handle.on_changed(update_plots)

# Setup for the amplitude slider of the second signal component
amp2_ax = plt.axes([0.08, 0.20, 0.22, 0.02])
amp2_handle = widgets.Slider(amp2_ax, 'Amplitude 2', 0, A_max, valinit=0)
amp2_handle.on_changed(update_plots)

# Setup for the phase slider of the second signal component
phase2_ax = plt.axes([0.08, 0.17, 0.22, 0.02])
phase2_handle = widgets.Slider(phase2_ax, 'Phase 2 [rad]', 0., 2*np.pi, valinit=0)
phase2_handle.on_changed(update_plots)

# Setup for the noise slider
noise_ax = plt.axes([0.08, 0.12, 0.22, 0.02])
noise_handle = widgets.Slider(noise_ax, 'Noise', 0, A_max, valinit=0)
noise_handle.on_changed(update_plots)

# Setup for the time windowing lower limit slider
t0_ax = plt.axes([0.08, 0.09, 0.22, 0.02])
t0_handle = widgets.Slider(t0_ax, 'T lower lim. [s]', 0, T, valinit=0)
t0_line, = signal_ax.plot([0,0],[-A_max,A_max],'k--')
t0_handle.on_changed(update_plots)

# Setup for the time windowing upper limit slider
t1_ax = plt.axes([0.08, 0.06, 0.22, 0.02])
t1_handle = widgets.Slider(t1_ax, 'T upper lim. [s]', 0, T, valinit=T)
t1_line, = signal_ax.plot([T,T],[-A_max,A_max],'k--')
t1_handle.on_changed(update_plots)

# Setup for the signal type radio buttons
waveform_ax = plt.axes([0.34, 0.16, 0.07, 0.1])
waveform_handle = widgets.RadioButtons(waveform_ax,('Sine','Square', 'Sawtooth'))
waveform_handle.on_clicked(update_plots)

# Setup for the frequency windowing lower limit slider
f0_ax = plt.axes([0.42, 0.42, 0.22, 0.02])
f0_handle = widgets.Slider(f0_ax, 'f lower limit [Hz]', 0, 100, valinit=0)
f0_line, = power_ax.plot([0,0],[0,np.max(power_spectrum)],'k--')
f0_handle.on_changed(update_plots)

# Setup for the frequency windowing upper limit slider
f1_ax = plt.axes([0.42, 0.39, 0.22, 0.02])
f1_handle = widgets.Slider(f1_ax, 'f upper limit [Hz]', 0, 100, valinit=50)
f1_line, = power_ax.plot([nyquist,nyquist],[0,np.max(power_spectrum)],'k--')
f1_handle.on_changed(update_plots)

plt.show()
