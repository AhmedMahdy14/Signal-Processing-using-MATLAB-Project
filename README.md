# Signal-Processing-using-MATLAB-Project
Part I: Sound Signal
1. Import the two songs into MATLAB and play them.
2. Increase the sample rate of the song "long track" to double its original value. Play the
sound and Save the new song (after changing the sample rate) to a .wav le (with a name
of your choice).
3. Reduce the sample rate of the song "short track" to one third of its original value. Play
the sound and Save the new song (after changing the sample rate) to a .wav le (with a
name of your choice).
4. Fade in the rst 2 seconds of the song "long track" and Fade out the last 2 seconds of the
same song. Play the sound.
5. Use a sine signal to switch between playing the 2 songs alternatively according to the sign
of the sine signal.
i.e
 Trim the length of the "long track" to be the same as the length of the "short track".
 Construct a sine signal of length similar to that of each of the two songs. Choose the
frequency to be 0.5 Hz.
 For the positive parts of the sine signal, play the song "long track".
 For the negative parts of the sine signal, play the song "short track".

Part II: General Signal Operations
2.2.1 It is required to implement a general signal generator that has the following specifications:
 When this part if the code runs, the program asks the user for the following parameters:
1. Sampling rate of the signal.
2. Start and end of time basis.
3. Number of the break points and their positions (i.e. the points that the signal defi
nitio rule changes).
Example: The signal is dened from -2 to 0 as a DC signal and from 0 to 2 as ramp
the user will enter that the number of break points is 1 and the position at t = 0.
 According to the number of break points the program asks the user at each region to enter
the specications of the signal at this region, which are:
1. For DC signal: Amplitude.
2. For Ramp signal: slope and intercept.
3. For Exponential signal: Amplitude and exponent.
4. For Sinusoidal signal: Amplitude, frequency and phase.
 Then the program will form the signal using the indirect method (review section 2 for the
indirect method).
 Plot the signal in time domain and frequency domain.


Create a system:
It is required to use the signal generated previously as an input to a system and get the response
(output) of the system due to this input.
How to create the system:
 First the program will ask the user to choose if he wants to dene the system via impulse
response or by a transfer function.
 If the user chooses to enter the impulse response:
1. The program should again form the impulse response by the same way of generation
of signal.
2. The program then calculates the output of the system (using the command "conv").
3. Plot the output in time domain and frequency domain.
 If the user chooses to enter the transfer function:
1. The user should enter the numerator and denominator of the transfer function of the
system.
2. The program should plot the z-plane of this system and specify whether it is stable,
unstable or marginally stable.
3. Then the program calculates the output of the system (using the command "filter" ).
4. Plot the output in time domain and frequency domain.
