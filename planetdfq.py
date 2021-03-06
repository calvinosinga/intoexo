import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve

class Planet():
    """
    The implementation of the model for the interior of exoplanets.
    """

    def __init__(self, materials = ['H_2'], mass_fractions = [1], 
            verbose = 0):
        """
        Create an instance of a planet interior model.
        
        TODO: currently the ability to model an interior with
        multiple materials with different mass fractions is not
        implemented, so stick to just one material for now.

        Args:
            materials (list, optional): The materials that are in 
                this planet's interior, in order from center to
                surface. Defaults to ['H_2'], for ideal gas testing 
                case.
            mass_fractions (list, optional): The fraction of the
                planet's mass that belongs to each material listed
                in materials. Defaults to [1].
            verbose (int, optional): If true, will print out
                in-progress statements during MR calculation loop.
                When testing, would be good to set to 1 (True).
                Defaults to 0.
        """
        if isinstance(materials, str):
            self.materials = [materials]
        else:
            self.materials = materials
        self.mfs = mass_fractions
        self.v = verbose
        return
    
    def eos(self, material, pressure, temp):
        """
        Calculates the density given a material, which is used to
        identify the equation of state to use, and a pressure and
        temperature.

        Args:
            material (string): the string id for the material to
                a particular eos. To implement new eos, use an
                if statement in this method!

            pressure (float): the pressure in N / m^2
            temp (float): the temperature in K

        Raises:
            NotImplementedError: if a given material does not
            yet have an eos, raises this error.

        Returns:
            rho: density in kg / m^3
        """

        if material == 'H_2':
            k_b = 1.38e-23 # J / K
            mol_weight = 2.016 # amu
            amu_to_kg = 1.66e-27 # kg / amu
            mol_weight *= amu_to_kg # now in kg
            return pressure * mol_weight/ k_b / temp

            
        elif material == 'H2O':
            rho0 = 1460
            c = 0.00311
            n = 0.513
            return rho0 + c*pressure**n
            
        elif material == 'MgSiO3': # perovskite
            rho0 = 4100
            c = 0.00161
            n = 0.541
            return rho0 + c*pressure**n
        
        elif material == 'Fe': # alpha
            rho0 = 8300
            c = 0.00349
            n = 0.528
            return rho0 + c*pressure**n
        
        else:
            msg = '%s not an implemented material'%material
            raise NotImplementedError(msg)
        
    
    def dpdr(self, r, mass, rho):
        """
        Calculate dp/dr

        Args:
            r (float): radius in m
            mass (float): mass of the shell in kg
            rho (float): density of the sphere in kg / m^3

        Returns:
            dpdr (float): dp / dr
        """
        G = 6.67e-11 # in SI


        return -G * mass * rho / r**2
    
    def dmdr(self, r, rho):
        """
        Calculate dm/dr.

        Args:
            r (float): radius in m
            rho (float): density in kg / m^3

        Returns:
            dmdr (float): dm / dr
        """

        return r**2 * 4 * np.pi * rho
    
    def shellMassFraction(self):
        # TODO: use to implement mass fraction later
        return 0

    def getMaterial(self, r):
        # TODO implement mass fraction
        return self.materials[0]

    def calculateMP(self, central_pressure, central_mass,
            rstart, rstop, temp, steps = int(1e6), stop_pressure = 0,
            step_type = 'log'):
        """ Calculates the mass and pressure as a function of radius.

        Args:
            central_pressure (float): the pressure at the center of the
                planet.
            central_mass (float): the mass at the center of the planet.
                Should typically be 0.
            rstart (float): Where to start the calculation from.
            rstop (float): a guess for where the radius will end.
                This will not determine when the computation ends,
                but rather serves as a way to take efficient dr steps.
            temp (func, float): If a radially varying temperature, pass in
                a function that only takes r as an input. Otherwise, can
                give a float and then a function will be created.
            steps (int, optional): max number of steps to take. Specify
                if likely that the integration will take too long.
                Defaults to int(1e6).
            stop_pressure (float, optional): The condition to stop the
                loop. When P reaches this value, computation stops.
                Defaults to 0.

        Returns:
            radii (array): the radii for the mass and pressure
                calculations.
            mass (array): the masses within the corresponding 
                radii.
            pressure (array): the pressures at the corresponding
                radii. 
        """
        # handle default behavior for temp is non-function
        if isinstance(temp, int) or isinstance(temp, float):

            def _oneVal(r):
                return temp
            
            temp_func = _oneVal
        else:
            temp_func = temp
        
        # function to handle ode integration
        def calcStep(pr_ma, r):
            P, m = pr_ma
            material = self.getMaterial(r)
            rho = self.eos(material, P, temp_func(r))
            dmdr = self.dmdr(r, rho)
            dpdr = self.dpdr(r, m, rho)
            return [dpdr, dmdr]

        
        idx = 0
        pressure = np.zeros(steps - 1)
        mass = np.zeros(steps - 1)
        radii = np.zeros(steps - 1)

        # numpy arrays run much, much faster than lists, 
        # but we don't know how long we need to loop before
        # the pressure condition is reached. So we create a numpy
        # array of arbitrary length and keep padding it until
        # finally the stopping pressure is reached.

        while pressure[-1] > stop_pressure or idx == 0:
            
            # the radii for this iteration
            if step_type == 'log':
                new_radii = np.geomspace(rstart, rstop, steps)
            elif step_type == 'lin':
                new_radii = np.linspace(rstart, rstop, steps)
            # if not the first iteration, need to pad RMP arrays
            # to store the new values
            if not idx == 0:
                pressure = np.pad(pressure, ((0, steps - 1)), constant_values = (0, 0))
                mass = np.pad(mass, ((0, steps - 1)), constant_values = (0, 0))
                radii = np.pad(radii, ((0, steps - 1)), constant_values = (0, 0))
            else:
                p0 = central_pressure
                m0 = central_mass

            out = odeint(calcStep, [p0, m0], new_radii)
            # since the last value is saved as initial value of
            # next iteration, exclude it.
            pressure[idx:] = out[:-1, 0]
            mass[idx:] = out[:-1, 1]
            radii[idx:] = new_radii[:-1]

            
            # get new r-range the same difference in logspace
            if step_type == 'log':
                dif = np.log10(new_radii[-1] / new_radii[0])
                tmp = rstop
                rstop = rstop * 10**(dif)
                rstart = tmp
            elif step_type == 'lin':
                tmp = rstop
                rstop *= 2
                rstart = tmp
            # set new initial params
            p0 = out[-1, 0]
            m0 = out[-1, 1]

            idx += steps - 1

        
        # there's a decent likelihood that we overshot the
        # stopping pressure, this trims excess values...
        stopidx = (np.abs(pressure - stop_pressure)).argmin()

        radii = radii[:stopidx]
        mass = mass[:stopidx]
        pressure = pressure[:stopidx]

        return radii, mass, pressure

    def integrateCP(self, central_pressures, central_mass,
            rstart, rstop, temp, steps = int(1e6), 
            stop_pressure = 0, step_type = 'log'):
        """
        Given an array of central pressures, will calculate the final
        mass and radius of the planet and return them for each of the
        given central pressures.

        The arguments are used as input into calculateMRP(...), see
        that function for more details
        Args:
            central_pressures (np.array): 1D array containing the
                central pressures (Pa). Used as input into
                calculateMRP(...).
            
            #### BELOW SEE CALCULATEMRP ####
            central_mass (float):
            rstart (float):
            rstop (float): 
            temp (func): 
            steps (int, optional):
            stop_pressure (int, optional): Defaults to 0.

        Returns:
            surface_radii: (np.array) 1D array containing
                the `surface` radius at which the pressure
                reached the stop condition for each of the
                given central pressures (m).
            
            total_masses: (np.array) 1D array that contains
                the mass enclosed in the surface radius (kg).
        """
        
        cps = central_pressures # for convenience
        surface_radii = np.zeros_like(cps)
        total_masses = np.zeros_like(cps)
        num = len(cps)
        # for each central pressure, calculate the final mass and
        # radius for when the stop condition is met
        for i in range(len(cps)):
            if i % int(num / 5) == 0:
                print('finished 1/5 of central pressures...')

            # pressure value not used
            r, m, _ = self.calculateMP(cps[i], central_mass, rstart, rstop, temp,
                    steps, stop_pressure, step_type)
            
            # save final radius and mass
            surface_radii[i] = r[-1]
            total_masses[i] = m[-1]
        
        return surface_radii, total_masses
    
    def calculateMPList(self, central_pressure, radii, temp, 
            max_steps = -1, stop_pressure = 0):
        """
        Calculates the mass, radius, and pressure for each spherical
        shell dr apart with the given temperature and central 
        pressure.  *OLD VERSION*

        Args:
            central_pressure (float): given central pressure value
                for P(r=0).

            radii (array): values to evaluate pressure and mass at.

            temp (function, int, float): specify the temperature.
                Can be given as a function that takes one input, radius,
                (i.e. write it like def temp(r): ... and pass it as a 
                parameter) to allow for a radially varying temperature.

            max_steps (int): sets a maximum number of steps, to prevent
                infinite loops
            
            stop_pressure (float): sets the pressure that serves as
                the end condition for the loop.

        Returns:
            masses (list): the mass of each spherical shell.

            pressures (list): the pressure at each spherical shell. 
        """

        # handle default behavior for temp is non-function
        if isinstance(temp, int) or isinstance(temp, float):

            def _oneVal(r):
                return temp
            
            temp_func = _oneVal
        else:
            temp_func = temp

        pressures = [central_pressure]
        masses = [0]
        radii = np.insert(radii, 0, 0)
        step = 1
        material_idx = 0
        
        # this loop integrates over a spherical shell of width dr
        # until the stop condition is met.
        # the condition to stop the loop as described in Seager is
        # when the pressure is 0 (2nd to last paragraph in S. 2)
        while pressures[-1] > stop_pressure and step < len(radii):
            if self.v:
                print('starting step %d'%step)
            elif step % 1000 == 0:
                print('starting step %d'%step)
            

            # TODO predict final mass and use it to determine the shell's
            # mass fraction, used to determine if we need to switch materials
            if self.shellMassFraction() > self.mfs[material_idx]:
                material_idx += 1
            shell_material = self.materials[material_idx]
            
            dr = radii[step] - radii[step - 1]
            
            # calculate rho using the previous pressure value using the material's
            # equation of state
            rho = self.eos(shell_material, pressures[-1], 
                    temp_func(radii[step]))

            if self.v:
                print('\tdensity value: ' + str(rho))

            # calculate the dm/dr, multiply by dr to get dm
            dm = self.dmdr(radii[step], rho) * dr

            if self.v:
                print('\tdm value calculated: ' + str(dm))
                print('\tshell mass: ' + str(masses[-1] + dm))


            # calculate dp/dr, multiply by dr to get dp
            
            dp = self.dpdr(radii[step], masses[-1], rho) * dr

            if self.v:
                print('\tdp value calculated: ' + str(dp))
                print('\tshell pressure: ' + str(pressures[-1] + dp))


            # save the values to be returned to respective lists
            masses.append(masses[-1] + dm)
            pressures.append(pressures[-1] + dp)

            # often dp is so small that this basically becomes
            # an infinite loop, so this is here to stop that
            # from breaking
            if step > max_steps and not max_steps == -1:
                print('warning: result not converged')
                break
                
            step += 1

        return np.array(masses), np.array(pressures)