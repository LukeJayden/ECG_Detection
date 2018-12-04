
# coding: utf-8

# # This is the ECG filtering 
# ## create a FIR filter 

# In[2]:


import numpy as np
import scipy
import pylab as plt

class FIR_filter:
    
    h=[]
    taps=0
    def __init__(self,_coefficientslist,taps):
        # code here
        for i in range(len(_coefficientslist)):
            self.h.append(_coefficientslist[i])
        self.taps=taps
        
    
    def dofilter(self,originalSignalSingle,nTh):
        # code here 
        i=0
        result=0
        while nTh-i>=0 and i< len(self.h):
            result+=self.h[i]*originalSignalSingle[nTh-i]
            i+=1
        return result
    
    @classmethod
    def bandstop_impulse_respouse_coefficients(cls):
        fs=1000
        f1 = 45.0/fs
        f2 = 55.0/fs
        n = np.arange(-200,200)
        h = (1/(n*np.pi))*(np.sin(f1*2*np.pi*n)-np.sin(f2*2*np.pi*n))
        h[200] = 1-(f1*2*np.pi-f2*2*np.pi)/np.pi;
        h = h * np.hamming(400)
        return h
    @classmethod
    def bandpass_impulse_respouse_coefficients(cls):
        fs=1000
        f1 = 45.0/fs
        f2 = 55.0/fs
        n = np.arange(-200,200)
        h = (1/(n*np.pi))*(np.sin(f2*2*np.pi*n)-np.sin(f1*2*np.pi*n))
        h[200] = (f1*2*np.pi-f2*2*np.pi)/np.pi;
        h = h * np.hamming(400)
        return h
    @classmethod
    def lowpass_impulse_respouse_coefficients(cls):
        fs=1000
        f1 = 45.0/fs
        f2 = 50.0/fs
        n = np.arange(-200,200)
        h = (1/(n*np.pi))*(np.sin(f2*2*np.pi*n))
        h[200] = f2*2;
        h = h * np.hamming(400)
        return h
    @classmethod
    def highpass_impulse_respouse_coefficients(cls):
        fs=1000
        
        f2 = 50.0/fs
        n = np.arange(-200,200)
        h = -(1/(n*np.pi))*(np.sin(f2*2*np.pi*n))
        h[200] = 1-f2*2;
        h = h * np.hamming(400)
        return h
    
    
# if __name__ == '__main__':
# #     h=[0.2,0.4,0.1,0.5]
#     data=np.loadtxt("2419351l.dat")
#     oringin=data[:,2]
# #     plt.plot(np.fft.fft(oringin))
#     h=FIR_filter.bandstop_impulse_respouse_coefficients()
# #     oringin=[1,2,1,1,0,0,0]
#     a=FIR_filter(h,len(h))
# #     plt.plot(oringin)
#     f=plt.figure()
#     plt.plot(a.h)
#     plt.title("bandstop")
#     f.savefig("bandstop.pdf",bbox_inches="tight")
    
# #     y=[]
# #     for i in range(len(oringin)):
# #         y.append(a.dofilter(oringin,i))
# #     plt.plot(y)
    
# #     print(a.h)
#     print(a.dofilter(oringin,1))


# In[3]:


import numpy as np
import scipy
import pylab as plt


class FIRfilter(FIR_filter):
    taps = 0
    Fstop = 0
    Fpass = 0
    Fs = 0
    h = []

    def __init__(self, Fs, choice, *args):
        if args[0] != None :
            self.Fstop = args[0]
        if len(args) == 2 :
            self.Fpass = args[1]
        self.Fs = Fs
        
        self.taps = int(Fs / (self.Fstop - self.Fpass) * 4)
        #         int(60/(22*((self.Fstop-self.Fpass)/self.Fs)))
        if choice == "bandpass":
            self.h = self.bandpass_impulse_respones(self.Fpass, self.Fstop, self.Fs)
        if choice == "bandstop":
            self.h = self.bandstop_impulse_respones(self.Fpass, self.Fstop, self.Fs)
        if choice == "lowpass":
            self.h = self.lowpass_impulse_respones(self.Fstop, self.Fs)
        if choice == "highpass":
            self.h = self.highpass_impulse_respones(self.Fstop, self.Fs)
        

    def bandstop_impulse_respones(self, freq1, freq2, Fs):
        f_resp = np.ones(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        f_resp[(int(freq1 / (Fs / self.taps)) - 1):(int(freq2 / (Fs / self.taps)) + 2)] = 0
        f_resp[int((self.taps - (freq2 / (Fs / self.taps)))) - 1:int((self.taps - (freq1 / (Fs / self.taps)))) + 2] = 0
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def bandpass_impulse_respones(self, freq1, freq2, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        f_resp[(int(freq1 / (Fs / self.taps)) - 1):(int(freq2 / (Fs / self.taps)) + 2)] = 1
        f_resp[int((self.taps - (freq2 / (Fs / self.taps)))) - 1:int((self.taps - (freq1 / (Fs / self.taps)))) + 2] = 1
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def lowpass_impulse_respones(self, frequency, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        f_resp[0:(int(frequency / (Fs / self.taps)) + 1)] = 1
        f_resp[
        int((self.taps - (frequency / (Fs / self.taps)))) - 1:int((self.taps - (frequency / (Fs / self.taps)))) + 2] = 1
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind

    def highpass_impulse_respones(self, frequency, Fs):
        f_resp = np.zeros(self.taps)
        # note we need to add "+1" to the end because the end of
        # the range is not included.
        f_resp[(int(frequency / (Fs / self.taps)) + 1):self.taps] = 1
        f_resp[
        int((self.taps - (frequency / (Fs / self.taps)))) - 1:int((self.taps - (frequency / (Fs / self.taps)))) + 2] = 1
        hc = np.fft.ifft(f_resp)
        h = np.real(hc)
        h_shift = []
        # this is from index 0 to index 49 on the left
        # and on the right hand side from index 50 to index 99
        h_shift[0:int(self.taps / 2)] = h[int(self.taps / 2):self.taps]
        h_shift[int(self.taps / 2):self.taps] = h[0:int(self.taps / 2)]
        h_wind = h_shift * np.hamming(self.taps)
        return h_wind



if __name__ == '__main__':
    a = FIRfilter(1000, "bandstop", 60, 45)
    data = np.loadtxt("ecg_2.dat")
    oringin = data[:, 2]
    oringin=oringin[15000:40000]
    y = []
    for i in range(len(oringin)):
        y.append(a.dofilter(oringin, i))
#     data[:, 2]=y
    f=plt.figure()
#     data[:, 2]
    plt.plot(oringin)
    plt.title("ECG")
    plt.show()
    plt.plot(y)
    plt.title("filtered ECG")
    plt.show()
#     plt.plot(oringin)
#     plt.show()

