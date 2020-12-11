**Description**

STEMMUS-UEB V1.0.0 coupled the Simultaneous Transfer of Energy, Momentum, and Mass in Unsaturated Soil with Freeze-Thaw (STEMMUS-FT, Zeng et al., 2011a,b; Zeng and Su, 2013; Yu et al., 2018) with snowmelt process model (Utah Energy Balance, UEB, Tarboton and Luce, 1996).

With the various representations of soil physical processes, from the basic coupled, to tightly advanced coupled water and heat transfer, and further to the explicit consideration of air flow, STEMMUS-FT model facilitates us to understand and interpret the role of soil physical processes (Zeng et al., 2011a,b; Yu et al., 2018; Yu et al., 2020).

UEB model is a physically based snowmelt model developed by David Tarboton's group. The model uses the lumped representation of snowpack with two primary states variables, snow water equivalent SWE and the internal energy U. Water and energy balance are numerically solved and kept track. UEB uses the physical based calculations of surface energy exchanges. Melt outflow is regarded as the function of the liquid fraction, by Darcy's law. The latest version, UEBGrid, developed for gridded application of UEB can be found from https://github.com/dtarb/UEBFortran. Older versions of UEB may be obtained from http://hydrology.usu.edu/dtarb/snow/snow.html.

The one-way sequential coupling is employed to couple the soil model (STEMMUS-FT) with the snowpack model (UEB). UEB model takes the atmospheric forcing as the input (precipitation, air temperature, wind speed and direction, relative humidity, shortwave and longwave radiation) and solves the snowpack energy and mass balance, provides the melt water flux and heat flux as the surface boundary conditions for the soil model STEMMUS-FT. STEMMUS-FT then solves the energy and mass balance equations of soil layers in one timestep.

**Setup and requirements**

The code is tested with MATLAB 2019b.
STEMMUS-UEB is executed in MATLAB by simply running `MainLoop.m` after you finish all the model setup and give the input data to STEMMUS-UEB. Several steps are necessary to build up the model setup.
-	Setting the temporal information and model domain;
-	Setting soil properties and snow properties;
-	Setting the initialization condition for soil and snow submodules, respectively;
-	Inputting the meteorological forcing information;
-	Setting the surface/bottom conditions;
Then you are ready to run STEMMUS-UEB by running `MainLoop.m`.
