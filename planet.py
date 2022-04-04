import numpy as np
from scipy.integrate import quad
import time

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
            mol_weight *= amu_to_kg
            return pressure * mol_weight/ R / temp
        
        else:
            msg = '%s not an implemented material'%material
            raise NotImplementedError(msg)
        
    
    def integratePressure(self, r1, r2, rho):
        """
        Calculate pressure within spherical shell

        Args:
            r1 (float): lower bounds of integral
            r2 (float): upper bounds of integral
            rho (float): density of the shell

        Returns:
            int_result (list): output of quad(...)
        """
        G = 6.67e-11 # in SI

        def _pressure(r):
            # Get m(r) by integrating (we do not assume constant mass
            # within a shell). m(r) is defined to be mass within r, 
            # so integrate from zero.
            mass = self.integrateMass(0, r, rho)

            # this is a version I tried without using the rho value, 
            # just the integrated mass value. Still didn't work.
            # V = 4/3 * np.pi * r**3
            # return mass[0]**2 / V * -G / r**2
            
            return mass[0] * rho * -G / r**2

        int_result = quad(_pressure, r1, r2)
        return int_result
    
    def integrateMass(self, r1, r2, rho):
        """
        Calculate the mass within a spherical shell.

        Args:
            r1 (float): lower bounds of integral
            r2 (float): upper bounds of integral
            rho (float): density of the shell

        Returns:
            int_result (list): output of quad(...)
        """
        def _mass(r):
            return 4 * np.pi * r**2 * rho
        
        int_result = quad(_mass, r1, r2)
        return int_result
    
    def shellMassFraction(self):
        # TODO: use to implement mass fraction later
        return 0

    def calculateMRP(self, central_pressure, dr, temp):
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

        # define values needed for the loop
        radii = [0]
        pressures = [central_pressure]
        masses = [0]
        step = 0
        material_idx = 0
        
        # this loop integrates over a spherical shell of width rstep
        # until the stop condition is met.
        # the condition to stop the loop as described in Seager is
        # when the pressure is 0 (2nd to last paragraph in S. 2)
        while pressures[-1] > 0:
            if step % 10 == 0 and self.v:
                print('starting step %d'%step)
            r1 = radii[-1]
            r2 = r1 + dr

            # TODO predict final mass and use it to determine the shell's
            # mass fraction, used to determine if we need to switch materials
            if self.shellMassFraction() > self.mfs[material_idx]:
                material_idx += 1
            

            # we assume that the dr is small enough that the density changes
            # slowly
            # we need this assumption since
            # P(r, m(r)) -> m(r, rho(r)) -> rho(r, P(r)) -> P(r, m(r)) -> ...
            # will cause just a circular integral
            # thus we use the last pressure value, which would be the value at
            # the boundary between this shell and the previous shell

            # calculate rho using the previous pressure value using the material's
            # equation of state
            rho = self.eos(self.materials[material_idx], 
                    pressures[-1], temp_func(r1))

            if self.v:
                print('\tdensity value: ' + str(rho))

            
            # calculate the mass by integrating over the conservation of mass eqtn
            if self.v:
                print('\tstarting mass integration...')

            start = time.time()
            mass = self.integrateMass(r1, r2, rho)

            if self.v:
                print('\tmass integration complete: %d s'%(time.time() - start))
                print('\tmass value calculated: ' + str(mass[0]))
                print('\tintegral abserr: ' + str(mass[1]))



            # calculate the pressure by integrating over the hydrostatic eqbm eqtn 
            if self.v:    
                print('\tstarting pressure integration...')
            
            start= time.time()
            pressure = self.integratePressure(r1, r2, rho)

            if self.v:
                print('\tpressure integration complete: %d s'%(time.time() - start))
                print('\tpressure value calculated: ' + str(pressure[0]))
                print('\tintegral abserr: ' + str(pressure[1]))



            # save the values to be returned to respective lists
            radii.append(r2)
            masses.append(mass[0] + np.sum(masses))
            pressures.append(pressure[0])

        return masses, radii, pressures

