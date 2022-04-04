import numpy as np

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
            R = 8.314 # N m / k / mol
            mol_weight = 2.016 # amu / mol
            amu_to_kg = 1.66e-27 # kg / amu
            mol_weight *= amu_to_kg # now in kg / mol
            return pressure * mol_weight/ R / temp
        
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
        Calculate the mass within a spherical shell.

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

    def calculateMRP(self, central_pressure, dr, temp, 
            max_steps = -1, stop_pressure = 0):
        """
        Calculates the mass, radius, and pressure for each spherical
        shell dr apart with the given temperature and central 
        pressure.
        Args:
            central_pressure (float): given central pressure value
                for P(r=0).

            dr (float): width of spherical shells.

            temp (function, int, float): specify the temperature.
                Can be given as a function that takes one input, radius,
                (i.e. write it like def temp(r): ... and pass it as a 
                parameter) to allow for a radially varying temperature.

        Returns:
            masses (list): the mass of each spherical shell.

            radii (list): the radius of each spherical shell.

            pressures (list): the pressure at each spherical shell. 
        """

        # handle default behavior for temp is non-function
        if isinstance(temp, int) or isinstance(temp, float):

            def _oneVal(r):
                return temp
            
            temp_func = _oneVal
        else:
            temp_func = temp

        # define initial values needed for the loop
        radii = [dr]
        pressures = [central_pressure]
        # initial mass value is m(r_1) = m(dr) = 0 + dm
        masses = [dr * self.dmdr(radii[0], self.eos(self.materials[0],
                central_pressure, temp_func(radii[0])))]
        step = 0
        material_idx = 0
        
        # this loop integrates over a spherical shell of width rstep
        # until the stop condition is met.
        # the condition to stop the loop as described in Seager is
        # when the pressure is 0 (2nd to last paragraph in S. 2)
        while pressures[-1] > stop_pressure:
            if self.v:
                print('starting step %d'%step)
            elif step % 1000:
                print('starting step %d'%step)
            
            r_new = radii[-1] + dr

            # TODO predict final mass and use it to determine the shell's
            # mass fraction, used to determine if we need to switch materials
            if self.shellMassFraction() > self.mfs[material_idx]:
                material_idx += 1
            

            # calculate rho using the previous pressure value using the material's
            # equation of state
            rho = self.eos(self.materials[material_idx], 
                    pressures[-1], temp_func(r_new))

            if self.v:
                print('\tdensity value: ' + str(rho))


            # calculate the dm/dr, multiply by dr to get dm
            dm = self.dmdr(r_new, rho) * dr
            if self.v:
                print('\tdm value calculated: ' + str(dm))
                print('\tshell mass: ' + str(masses[-1] + dm))


            # calculate dp/dr, multiply by dr to get dp
            
            dp = self.dpdr(r_new, masses[-1] + dm, rho) * dr

            if self.v:
                print('\tdp value calculated: ' + str(dp))
                print('\tshell pressure: ' + str(pressures[-1] + dp))


            # save the values to be returned to respective lists
            radii.append(r_new)
            masses.append(masses[-1] + dm)
            pressures.append(pressures[-1] + dp)

            if step > max_steps and not max_steps == -1:
                break
                
            step += 1

        return masses, radii, pressures

