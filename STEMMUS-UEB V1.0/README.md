STEMMUS-UEB V1.0 coupled the Simultaneous Transfer of Energy, Momentum, and Mass in Unsaturated Soil (STEMMUS) with snowmelt process model (Utah Energy Balance, UEB). 

With the various representations of soil physical processes, from uncoupled, to tightly coupled water and heat transfer, and further to the explicit consideration of air flow, 
STEMMUS model facilitate us to understand and interpret the role of soil physcial processes. 
In the curret version, STEMMUS further takes into account the freeze/thaw (FT) process and FT induced water and heat coupling transfer (STEMMUS-FT). 

UEB model is physcially based snowmelt model developed by David Tarboton's group. The model uses the lumped representation of snowpack with two primary states variables,
snow water equivalent W and the internal energy U. Water and energy balance are numerically solved and kept track. UEB uses the physical based calculations of surface energy exchanges. Melt outflow is regarded as the function of the liquid fraction, by Darcy's law. The latest version, UEBGrid, developed for gridded application of UEB can be found from https://github.com/dtarb/UEBFortran. Older versions of UEB may be obtained from http://hydrology.usu.edu/dtarb/snow/snow.html

The one-way sequential coupling is employed to couple the soil model (STEMMUS) with the snowpack model (UEB). UEB model takes the atmospheric forcing as the input (precipitation, air temperature, wind speed and direction, relative humidity, shortwave and longwave radiation) and solves the snowpack energy and mass balance, provides the melt water flux and heat flux as the surface boundary conditions for the soil model STEMMUS. STEMMUS then solves the energy and mass balance equations of soil layers in one timestep.

