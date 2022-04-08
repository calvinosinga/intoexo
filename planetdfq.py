import numpy as np
from scipy.integrate import odeint

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
        if material == 'H_2': # using ideal gas law
            k_b = 1.38e-23 # J / K
            mol_weight = 2.016 # amu
            amu_to_kg = 1.66e-27 # kg / amu
            mol_weight *= amu_to_kg # now in kg
            return pressure * mol_weight/ k_b / temp
        
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
            rstart, rstop, temp, steps = int(1e6), stop_pressure = 0):

        # handle default behavior for temp is non-function
        if isinstance(temp, int) or isinstance(temp, float):

            def _oneVal(r):
                return temp
            
            temp_func = _oneVal
        else:
            temp_func = temp
        
        def calcStep(pr_ma, r):
            P, m = pr_ma
            material = self.getMaterial(r)
            rho = self.eos(material, P, temp_func(r))
            dmdr = self.dmdr(r, rho)
            dpdr = self.dpdr(r, m, rho)
            return [dpdr, dmdr]

        
        idx = 0
        pressure = np.zeros(steps)
        mass = np.zeros(steps)
        radii = np.zeros(steps)

        pressure[0] = central_pressure
        mass[0] = central_mass

        while pressure[-1] > stop_pressure or idx == 0:
            

            new_radii = np.geomspace(rstart, rstop, steps)

            if not idx == 0:
                pressure = np.pad(pressure, ((0, steps)), constant_values = (0, 0))
                mass = np.pad(mass, ((0, steps)), constant_values = (0, 0))
                radii = np.pad(radii, ((0, steps)), constant_values = (0, 0))
                
            out = odeint(calcStep, [pressure[idx], mass[idx]], new_radii)
            
            
            pressure[idx:] = out[:, 0]
            mass[idx:] = out[:, 1]
            radii[idx:] = new_radii[:]

            temp = rstop
            rstop = (np.log10(rstop) - np.log10(rstart)) * steps
            rstart = temp

            idx += steps

        mass = np.trim_zeros(mass, 'b')
        pressure = np.trim_zeros(pressure, 'b')
        radii = radii[:pressure.shape[0]]
        return radii, mass, pressure

    def integrateCP(self, central_pressures, central_mass,
            rstart, rstop, temp, stop_pressure = 0):
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
            dr (float): 
            temp (func): 
            max_steps (int, optional): Defaults to -1.
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
        
        # for each central pressure, calculate the finall mass and
        # radius for when the stop condition is met
        for i in range(len(cps)):

            # pressure value not used
            r, m, _ = self.calculateMP(cps[i], central_mass, rstart, rstop, temp,
                    stop_pressure = stop_pressure)
            
            # save final radius and mass
            surface_radii[i] = r[-1]
            total_masses[i] = m[-1]
        
        return surface_radii, total_masses
    
    def calculateMPList(self, central_pressure, radii, temp, 
            max_steps = -1, stop_pressure = 0):
        """
        Calculates the mass, radius, and pressure for each spherical
        shell dr apart with the given temperature and central 
        pressure.

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