
# coding: utf-8

# # Heartrate Detection using FIR filter 
# 
# ## to detect momentary heart rate by measuring intervals of two heart beat

# In[35]:


import numpy as np
import scipy
import pylab as plt


class FIR_filter:
    h = []
    taps = 0

    def __init__(self, _coefficientslist, taps):
        # code here
        for i in range(len(_coefficientslist)):
            self.h.append(_coefficientslist[i])
        self.taps = taps

    def dofilter(self, signal):
        # code here
        i = 0
        result = 0
        #         result=0
        #         while nTh-i>=0 and i< len(self.h):
        #             result+=self.h[i]*originalSignalSingle[nTh-i]
        #             i+=1
        if self.taps !=0:
            for i in range(self.taps):
                result += self.h[i] * signal[i]
        else:
            for i in range(len(self.h)):
                result += self.h[i] * signal[i]


        return result

    @classmethod
    def bandstop_impulse_respouse_coefficients(cls):
        #sampling rate is 1000hz
        fs=1000
        #set up cutoff frequency
        f1 = 45.0/fs
        f2 = 55.0/fs
        n = np.arange(-200,300)
        #calculate the impulse response
        h = (1/(n*np.pi))*(np.sin(f1*2*np.pi*n)-np.sin(f2*2*np.pi*n))
        h[200] = 1-(f1*2*np.pi-f2*2*np.pi)/np.pi;
        #add the window function
        h = h * np.hamming(500)
        return h

    @classmethod
    def bandpass_impulse_respouse_coefficients(cls):
        #sampling rate is 1000hz
        fs=1000
        #set up cutoff frequency
        f1 = 45.0/fs
        f2 = 55.0/fs
        n = np.arange(-200,200)
        #calculate the impulse response
        h = (1/(n*np.pi))*(np.sin(f2*2*np.pi*n)-np.sin(f1*2*np.pi*n))
        h[200] = (f1*2*np.pi-f2*2*np.pi)/np.pi;
        #add the window function
        h = h * np.hamming(400)
        return h

    @classmethod
    def lowpass_impulse_respouse_coefficients(cls):
        #sampling rate is 1000hz
        fs=1000
        #set up cutoff frequency
        f1 = 45.0/fs
        f2 = 50.0/fs
        n = np.arange(-200,200)
        #calculate the impulse response
        h = (1/(n*np.pi))*(np.sin(f2*2*np.pi*n))
        h[200] = f2*2;
        #add the window function
        h = h * np.hamming(400)
        return h
    @classmethod
    def highpass_impulse_respouse_coefficients(cls):
        #sampling rate is 1000hz
        fs=1000
        #set up cutoff frequency
        f2 = 50.0/fs
        n = np.arange(-200,200)
        #calculate the impulse response
        h = -(1/(n*np.pi))*(np.sin(f2*2*np.pi*n))
        h[200] = 1-f2*2;
        #add the window function
        h = h * np.hamming(400)
        return h

# if __name__ == '__main__':
# #     h=[0.2,0.4,0.1,0.5]
#     data=np.loadtxt("2419351l.dat")
#     oringin=data[:,2]
# #     plt.plot(np.fft.fft(oringin))
#     h=FIR_filter.lowpass_impulse_respouse_coefficients()
# #     oringin=[1,2,1,1,0,0,0]
#     a=FIR_filter(h,len(h))
# #     plt.plot(oringin)
# #     plt.plot(a.h)

#     y=[]
#     for i in range(len(oringin)):
#         y.append(a.dofilter(oringin,i))
#     plt.plot(y)

#     print(a.h)
#     print(a.dofilter(oringin,1))


# In[36]:


import numpy as np
import scipy
import pylab as plt


class HeartrateDetection(FIR_filter):
    taps = 0
    Fstop = 0
    Fpass = 0
    Fs = 0
    h = []

    def __init__(self, Fs, choice, *args):
        if args[0] != None and choice != "detection":
            self.Fstop = args[0]
        if len(args) == 2 and choice != "detection":
            self.Fpass = args[1]
        self.Fs = Fs
        if choice != "detection":
            self.taps = int(Fs / (self.Fstop - self.Fpass) * 3)
        #         int(60/(22*((self.Fstop-self.Fpass)/self.Fs)))
        if choice == "bandpass":
            self.h = self.bandpass_impulse_respones(self.Fpass, self.Fstop, self.Fs)
        if choice == "bandstop":
            self.h = self.bandstop_impulse_respones(self.Fpass, self.Fstop, self.Fs)
        if choice == "lowpass":
            self.h = self.lowpass_impulse_respones(self.Fstop, self.Fs)
        if choice == "highpass":
            self.h = self.highpass_impulse_respones(self.Fstop, self.Fs)
        if choice == "detection":
            self.h = args[0]

    def bandstop_impulse_respones(self, freq1, freq2, Fs):
        '''
        function for generating bandstop impulse respones

        input:
                freq1: the start index of frequency
                freq2: the end index of frequency
                Fs: sampleing rate
        output:
                the coefficients of impulse response


        '''
        f_resp = np.ones(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        # Set the ideal frequency response of bandstop filter
        f_resp[(int(freq1 / (Fs / self.taps)) - 1):(int(freq2 / (Fs / self.taps)) + 2)] = 0
        f_resp[int((self.taps - (freq2 / (Fs / self.taps)))) - 1:int((self.taps - (freq1 / (Fs / self.taps)))) + 2] = 0
        # Inverse Fourier Transform
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        # shift it into pos.time
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        # Add the window function
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def bandpass_impulse_respones(self, freq1, freq2, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        # Set the ideal frequency response of bandpass filter
        f_resp[(int(freq1 / (Fs / self.taps)) - 1):(int(freq2 / (Fs / self.taps)) + 2)] = 1
        f_resp[int((self.taps - (freq2 / (Fs / self.taps)))) - 1:int((self.taps - (freq1 / (Fs / self.taps)))) + 2] = 1
        # Inverse Fourier Transform
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        # shift it into pos.time
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        # Add the window function
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def lowpass_impulse_respones(self, frequency, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        # Set the ideal frequency response of lowpass filter
        f_resp[0:(int(frequency / (Fs / self.taps)) + 1)] = 1
        f_resp[int((self.taps - (frequency / (Fs / self.taps)))) - 1:int((self.taps - (frequency / (Fs / self.taps)))) + 2] = 1
        # Inverse Fourier Transform
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        # shift it into pos.time
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        # Add the window function
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def highpass_impulse_respones(self, frequency, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        # Set the ideal frequency response of highpass filter
        f_resp[(int(frequency / (Fs / self.taps)) + 1):self.taps] = 1
        f_resp[int((self.taps - (frequency / (Fs / self.taps)))) - 1:int((self.taps - (frequency / (Fs / self.taps)))) + 2] = 1
        # Inverse Fourier Transform
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        # shift it into pos.time
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        # Add the window function
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    @classmethod
    def fir_coeff_generate(self, instance, TextUrl):
        '''

        '''
        data = np.loadtxt(TextUrl)
        oringin = data[:, 2]
        oringin = oringin[15000:40000]
        filteredSignal = instance.filtering(oringin)
        template = filteredSignal[16000:16600]
        #reverse the template 
        fir_coeff = template[::-1]
        return fir_coeff

    def filtering(self, oringin):
        y = []
        if self.taps != 0:
            signal = [0] * self.taps
        else:
            signal= [0] * len(self.h)
        for elem in oringin:
            signal.insert(0, elem)
            signal.pop(len(signal) - 1)

            y.append(self.dofilter(signal))
        #         for i in self.taps:
        #             y.append(self.dofilter())

        return y

    @classmethod
    def HeartrateDetect(cls, filteredSignal):
        start = []
        flag = 0
        for elem in filteredSignal:
            if elem > 20:
                start.append(list(filteredSignal).index(elem))
                continue
            if len(start) != 0:
                temp = np.median(start)
                start = []
                heartrate = 60 / (((temp - flag) / 25000) * 25)
                flag = temp
                print(int(heartrate))


if __name__ == '__main__':
    a = HeartrateDetection(1000, "bandstop", 60, 45)
    data = np.loadtxt("ecg_1.dat")
    oringin = data[:, 2]
    oringin = oringin[15000:40000]
    filteredSignal = oringin
    template = filteredSignal[16000:16600]
    fir_coeff = list(template[::-1])
    # fir_coeff = HeartrateDetection.fir_coeff_generate(a, "ecg_2.dat")
    # print(fir_coeff)
    data = np.loadtxt("ecg_2.dat")
    oringin = data[:, 3]
    oringin = oringin[15000:40000]
    b = HeartrateDetection(1000, "detection", fir_coeff)

    detorin = np.power(b.filtering(oringin),2)
    #     plt.plot(oringin)

    start = []
    flag = 0
    result = []
    for elem in detorin:
        if elem > 1.6:
            start.append(list(detorin).index(elem))
            continue
        if len(start) != 0:
            temp = np.median(start)
            start = []
            heartrate = 60 / (((temp - flag) / 25000) * 25)
            flag = temp
            result.append(int(heartrate))
    # plt.plot(detorin)
    plt.plot(result)
    plt.title("heartrate")
    plt.xlabel("time(s)")
    plt.show()

