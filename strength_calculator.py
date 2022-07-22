# Backfill Strength Calculator
"""
This file calculates the Backfill strength
based on Mitchell(1991) and 'Li & Aubertin (2012)' equations.

"""
import numpy as np
from colorama import Fore, Style

"""
Defining the stope geometry.
"""
### Exposed Face A1 ###
#Vertical = L
L = 16.8
#Horizontal length = H
H = 20
#Depth = B
B = 12.5

#Face Area (m^2) = f_area
f_area = 295.4
#Face Perimeter(m) = f_perimeter
f_perimeter = 70.2

# F/W Dip (degree) = fw
fw = 70
#HyD Radius (m)= Hyd_R 
Hyd_R = f_area/f_perimeter
#Fill lifts = fill_Lift
fill_Lift = 91.7
#Depth Incr (m) = D_incr
D_incr = 0.2
#Sill Pour Height (m)= H_sill_pour 
H_sill_pour = 6

'''
Geotechnical Parameters
'''
#alpha (degree) = Friction angle 
alpha = 34
#phi (degree)= Friction angle  (same as alpha?)
phi = 34
#Gravity (m/s^2)= g 
g = 9.81
#He = equivalent height of wedgeblock
He = H - (B* np.tan(alpha*np.pi/180))/2
#Fill Density (kg/m^3) = density
density = 2100
#Poisson ratio = poisson
poisson = 0.36
#factor of safety = fs
fs = 1.5
#tan(alpha) = tan_alpha
tan_alpha = np.tan(np.pi/180 *(45+alpha))

'''
Stress Profiles
'''
#unit weight = gama
gama = g * density/1000
#tan(phi) = tan_phi
tan_phi = np.tan(np.pi/180*phi)
#coeffK = K
K = 1/(1+2*(tan_phi**2))
#2nd HR = hr_2
hr_2 = H*B/(2*(H+B))
#sin(2*alpha) = sin_2alpha
sin_alpha = np.sin((45+alpha/2)*np.pi/180)
#sin(2*alpha) = sin_2alpha
sin_2alpha = np.sin(2*(45+alpha/2)*np.pi/180)

'''
Error Checks
'''
#H/B check
def hb_check():
    hb = H/B
    if hb < 5:
        hb_ch = 'OK!'
    else:
        hb_ch ='Error!'
    return hb_ch

print('H/B check is', Fore.CYAN + hb_check())
print(Style.RESET_ALL)

#Phi check
def phi_check():
    if phi < 40:
        p_ch = 'OK!'
    else:
        p_ch = 'Error!'
    return p_ch

print('Phi check is', Fore.CYAN +phi_check())
print(Style.RESET_ALL)

#LAR vs. HAR
def LAR_HAR_check():
    hb = H/B
    if hb > tan_alpha:
        lar_har = 'HAR'#High aspect ratio

    else:
        lar_har = 'LAR' #Low aspect ratio
    return lar_har
print('LAR vs. HAR check shows the stope is', Fore.CYAN + LAR_HAR_check())
print(Style.RESET_ALL)

#He check
def he_check():
    if He > 0:
        h_ch = 'OK!'
    else:
        h_ch ='Error!'
    return h_ch
print('He check is', Fore.CYAN + he_check())
print(Style.RESET_ALL)

'''
Paste Strength - Free Face 
'''
#self weight stress(kPa)= w_stress
w_stress = gama*H
#average sigma(v)
#FOS = fos
fos = 1

#LAR & HAR cohesions
def cohision():
    c_har = (gama* He/2)/((((fs - tan_phi/tan_alpha)*sin_2alpha)**-1)+He/L)
    c_lar = (gama* H/2)/(2*(((fs - tan_phi/tan_alpha)*sin_2alpha)**-1)+H/L)
    return c_har, c_lar
print('LAR cohesion =', round(cohision()[1],1), '\n HAR cohesion =', round(cohision()[0],1) )
print(Style.RESET_ALL)

#Simple Mitchell UCS
def mitchell():
    ucs_mitchell = gama* H/(1+H/L)
    return ucs_mitchell
print('Simple Mitchell UCS =',  round(mitchell(),1), 'kPa')
print(Style.RESET_ALL)

#Li & Auberton UCS (assumption: p0 = 0, rb = 1 )
def li_aubertin():
    for lar_har in LAR_HAR_check (): 
        if lar_har == 'HAR':
            c = (gama* He/2)/((((fs - tan_phi/tan_alpha)*sin_2alpha)**-1)+He/L)
        else:
            c = (gama* H/2)/(2*(((fs - tan_phi/tan_alpha)*sin_2alpha)**-1)+H/L)

    ucs_li_aubertin = 2*c*np.tan((45+phi/2)*np.pi/180)
    return ucs_li_aubertin
print('Li&Aubertin UCS =', round(li_aubertin(),1), 'kPa')
print(Style.RESET_ALL)

#Frictionless formula
def frictionless():
    f_less = fs* ((gama*L - cohision()[1])*(H-L/2)*tan_alpha)*sin_alpha/L
    return f_less
print('Frictionless UCS =', round(frictionless(),1), 'kPa')
print(Style.RESET_ALL)

##Design Strength
def design():
    expo_1 = mitchell() * fs
    expo_2 = frictionless()
    expo_3 = fs * w_stress
    return expo_1, expo_2, expo_3
print('Design Strengths for Exposure number=1 :', round(design()[0],1), 'kPa\n',
      'Design Strengths for Exposure number=2 :', round(design()[1]),1, 'kPa\n',
      'Design Strengths for Exposure number=3 :', round(design()[2],1), 'kPa')
