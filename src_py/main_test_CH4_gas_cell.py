# 1d interpolation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate as integr
import matplotlib.pyplot as plt
#
#------------------------------------------------------------------------------
fname_filter = 'filter.dat'
nu_min_filter = 0.4081632653E+04
dnu_filter = 0.3677146534E+00
npt_filter = 1155
nu_filter = np.linspace(nu_min_filter,
                        nu_min_filter+(npt_filter-1)*dnu_filter,
                        num=npt_filter)
data_filter = np.loadtxt(fname_filter, comments='#')
filter_crs =  np.ravel(data_filter)
finterp = interp1d(nu_filter, filter_crs, kind='cubic')#kind='cubic'
#------------------------------------------------------------------------------
fname_convol = 'convol.dat'
nu_min_convol = 0.4081901000E+04
dnu_convol = 0.002
npt_convol = 211900 
nu_convol = np.linspace(nu_min_convol,
                        nu_min_convol+(npt_convol-1)*dnu_convol,
                        num=npt_convol)
filter_convol = finterp(nu_convol)

data_convol = np.loadtxt(fname_convol, comments='#')
dagr_convol = np.ravel(data_convol)
#------------------------------------------------------------------------------
fname_ch4 = 'tau_abs_gcell_8cm_296K_1000mbar.txt'
data_ch4 = np.loadtxt(fname_ch4, comments='#')
tau = data_ch4[:, 1]
Trans = np.exp(-tau)
my_convol = filter_convol*Trans
#------------------------------------------------------------------------------
dagr_signal = integr.simps(dagr_convol, dx=dnu_convol)
my_signal = integr.simps(my_convol, dx=dnu_convol)
print('SIMPS: dagr signal = %.2f, my signal=%.2f, ratio=%.4f'%(dagr_signal, my_signal, dagr_signal/my_signal))
dagr_signal = integr.trapz(dagr_convol, dx=dnu_convol)
my_signal = integr.trapz(my_convol, dx=dnu_convol)
print('TRAPZ: dagr signal = %.2f, my signal=%.2f, ratio=%.4f'%(dagr_signal, my_signal, dagr_signal/my_signal))
#------------------------------------------------------------------------------
fname_S0 = '00_S0chkur.db'
data_S0 = np.loadtxt(fname_S0, skiprows=1)
nu_S0 = data_S0[:,0]
S0 = data_S0[:,1]
finterp_S0 = interp1d(nu_S0, S0, kind='cubic')#kind='cubic'
S0_ch4 = finterp_S0(nu_convol)  
#------------------------------------------------------------------------------
#
#fig.suptitle('Signal ratio: $\perp$ / $\parallel$', fontsize=12)
fig, axs = plt.subplots(nrows=4, ncols=1)
#fig.set_size_inches((12, 6))
#
axs[0].plot(nu_filter, filter_crs, 'b',
            nu_convol, filter_convol, ':r',
            nu_convol, dagr_convol, 'k',
            nu_convol, dagr_convol - my_convol, 'r')
axs[0].grid(True)
axs[0].set_ylabel('a.u')
axs[0].set_title('Signal', size=12)
axs[0].legend(['filter crs grid', 'filter fine grid', 'dagr convol', 'dagr-myself'])

axs[1].plot(nu_convol, tau, 'b')
axs[1].grid(True)
axs[1].set_xlabel('Wavenumber [cm-1]', fontsize=10)
axs[1].set_ylabel('$\\tau$')

axs[2].plot(nu_convol, Trans, 'g')
axs[2].grid(True)
axs[2].set_ylabel('exp(-$\\tau$)')

axs[3].plot(nu_convol, S0_ch4, 'm')
axs[3].grid(True)
axs[3].set_ylabel('S0($\\nu$), W/cm-1/cm2')
axs[3].set_xlabel('$\\nu$ [cm-1]', fontsize=10)