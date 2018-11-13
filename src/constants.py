""" model grid constants """

imt, jmt, km = 3600, 2400, 42     # POP high resolution grid dimensions


""" physical constants (SI units) """

g                 = 9.81          # [m s^-2]         gravitional constant
rho_sw            = 1.026         # [g cm^-3]        reference density of sea water
cp_sw             = 3996.         # [J kg^-1 K^-1]   heat capacity of sea water
R_earth           = 6.371e6       # [m]              radius of the Earth
A_earth           = 5.101e14      # [m^2]            surface area of the Earth
abs_zero          = -273.15       # [K]              absolute freezing
latent_heat_vapor = 2501000.      # [J kg^-1]        latent heat of vaporization


""" unit conversions """

spy               = 3600*24*365   # [s yr^-1]
Jpy_to_Wpsm       = 1/spy/A_earth # [J/yr] to [W/m^2]
