The KITE package for pre-processing is split up in various subclasses and contains several functions:

* [*class* `#!python kite.StructuralDisorder(lattice, concentration=0, position=None)` - Add disorder to the lattice.][structural_disorder]
* [*class* `#!python kite.Disorder(lattice)` - Add Guassian disorder to the lattice.][disorder]
* [*class* `#!python kite.Modification(magnetic_field=None, flux=None)` - Add a magnetic field to the lattice.][modification]
* [*class* `#!python kite.Configuration([...])` - Define the basic parameters used in the calculation.][configuration]
* [*class* `#!python kite.Calculation(configuration=None)` - Describe the required target functions.][calculation]
* [*function make_pybinding_model*][make_pybinding_model]
* [*function estimate_bounds*][estimate_bounds]
* [*function config_system*][config_system]
* [*warning LoudDeprecationWarning*][loud_deprecation_warning]

## StructuralDisorder

!!! declaration-class "*class* `#!python kite.StructuralDisorder(lattice, concentration=0, position=None)`"
    
     
:   Class that introduces Structural Disorder into the initially built [`#!python pb.Lattice`][lattice].

:   **Parameters**
    : <span id='structuraldisorder-lattice'>`#!python lattice`: *[`#!python pb.Lattice`][lattice]*</span>
        : The lattice used to build the structural disorder.
    : <span id='structuraldisorder-concentration'>`#!python concentration`: *`#!python float`*</span>
        : Concentration of disorder *(can only be different from `#!python 0` if `#!python position=None`)*.
    : <span id='structuraldisorder-position'>`#!python position`: *`#!python array_like`*</span>
        : Exact position of disorder *(can only be different from `#!python None` if `#!python concentration=0`)*.

:   **Methods**
    :   | Method                                                                                                         | Description                                                          |
        |----------------------------------------------------------------------------------------------------------------| -------------------------------------------------------------------- |
        | [`#!python add_vacancy(*disorder)`][structuraldisorder-add_vacancy]                                            | Add vacancy disorder to the lattice.                                 |
        | [`#!python add_structural_disorder(*disorder)`][structuraldisorder-add_structural_disorder]                    | Add structural disorder to the lattice.                              |
        | [`#!python add_local_vacancy_disorder(relative_index, sub)`][structuraldisorder-add_local_vacancy_disorder]    | Internal function to add one vacancy disorder to chosen position.    |
        | [`#!python add_local_bond_disorder(relative_index_from, [, ...])`][structuraldisorder-add_local_bond_disorder] | Internal function to add one bond disorder between chosen positions. |
        | [`#!python add_local_onsite_disorder(relative_index, [, ...])`][structuraldisorder-add_local_onsite_disorder]  | Internal function to add one onsite disorder to chosen position.     |
        | [`#!python map_the_orbital(orb, nodes_map)`][structuraldisorder-map_the_orbital]                               | Internal function to map the orbitals to the NodesMap.               |

    :   !!! declaration-function "<span id='structuraldisorder-add_vacancy'>*function* `#!python add_vacancy(*disorder)`"
            
            
        :   Add vacancy disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                  | Description                                                                                                       |
                |----------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
                | `#!python *disorder`:*`#!python str` or `#!python tuple(array_like, str)`* | Vacancy disorder, in the form of *`#!python sublattice_name`* or *`#!python ([relatice_index], sublattice_name)`* |
    
    :   !!! declaration-function "<span id='structuraldisorder-add_structural_disorder'>*function*`#!python add_structural_disorder(*disorder)`</span>"
            
            
        :   Add structural disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                                                                                                                                                                                       | Description                                                                                                                                                                                                            |
                |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python *disorder`:*`#!python tuple(array_like, str, float)` or `#!python tuple(array_like, str, array_like, str, array_like)` or `#!python tuple(array_like, str, array_like, str, float)` or `#!python tuple(array_like, str, array_like)`* | Vacancy disorder, in the form of *`#!python ([relatice_index], sublattice_name, onsite_energy)`* or *`#!python ([relatice_index_from], sublattice_name_from, [relatice_index_to], sublattice_name_to, onsite_energy)`* |
    
    

    :   !!! declaration-function "<span id='structuraldisorder-add_local_vacancy_disorder'>*function*`#!python add_local_vacancy_disorder(relative_index, sub)`</span>"
            
            
        :   Internal function to add one vacancy disorder to chosen position.
            
             **Parameters**

             :  | Parameter                                         | Description                                                 |
                |---------------------------------------------------|-------------------------------------------------------------|
                | `#!python relative_index`:*`#!python array_like`* | Relative index of the position to change the onsite energy. |
                | `#!python sub`:*`#!python str`*                   | Name of the sublattice to change the onsite energy.         |

    :   !!! declaration-function "<span id='structuraldisorder-add_local_bond_disorder'>*function* `#!python add_local_bond_disorder(relative_index_from, from_sub, relative_index_to, to_sub, hoppings)`</span>"
            
            
        :   Internal function to add one bond disorder between chosen positions.
            
             **Parameters**

             :  | Parameter                                                       | Description                                                                                                                        |
                |-----------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python relative_index_from`:*`#!python array_like`*          | Relative index of the position from wich the bond connects to change the onsite energy.                                            |
                | `#!python from_sub`:*`#!python str`*                            | Name of the sublattice from wich the bond connects to change the onsite energy.                                                    |
                | `#!python relative_index_to`:*`#!python array_like`*            | Relative index of the position to wich the bond connects to change the onsite energy.                                              |
                | `#!python to_sub`:*`#!python str`*                              | Name of the sublattice to wich the bond connects to change the onsite energy.                                                      |
                | `#!python hoppings`:*`#!python float` or `#!python array_like`* | The hopping energy between the different sublattices at the given positions, with the right shape to connect between the orbitals. |
    
    :   !!! declaration-function "<span id='structuraldisorder-add_local_onsite_disorder'>*function* `#!python add_local_onsite_disorder(relative_index, sub, value)`</span>"
            
            
        :   AInternal function to add one onsite disorder to chosen position.
            
             **Parameters**

            :  | Parameter                                                     | Description                                                                                    |
               |--------------------------------------------------------------|------------------------------------------------------------------------------------------------|
               | `#!python relative_index`:*`#!python array_like`*            | Relative index of the position to change the onsite energy.                                    |
               | `#!python sub`:*`#!python str`*                              | Name of the sublattice to change the onsite energy.                                            |
               | `#!python value`:*`#!python float` or `#!python array_like`* | The onsite energy of sublattice at the given positions, with the right shape for the orbitals. |
    
    :   !!! declaration-function "<span id='structuraldisorder-map_the_orbital'>*function* `#!python map_the_orbital(orb, nodes_map)`</span>"
            
            
        :   Internal function to map the orbitals to the NodesMap.
            
             **Parameters**

            :   | Parameter                              | Description                                                   |
                |----------------------------------------|---------------------------------------------------------------|
                | `#!python orb`:*`#!python str`*        | Name of the sublattice to give a unique value.                |
                | `#!python nodes_map`:*`#!python dict`* | The object that stores the unique values for the sublattices. |

## Disorder

!!! declaration-class "*class* `#!python kite.Disorder(lattice)`"
    
     
:   Class that introduces Disorder into the initially built [`#!python pb.Lattice`][lattice].

    The informations about the disorder are the *type*, *mean value*, and *standard deviation*.
    The function that you can use in the bulding of the [`#!python pb.Lattice`][lattice] is `#!python add_disorder()`.
    The class method takes care of the shape of the disorder chosen
    (it needs to be same as the number of orbitals at a given atom),
    and takes care of the conversion to the c++ orbital-only format.

:   **Parameters**
    : <span id='disorder-lattice'>`#!python lattice`: *[`#!python pb.Lattice`][lattice]*</span>
        : The lattice used to build the disorder.

:   **Methods**
    :   | Method                                                                                | Description                                              |
        | ------------------------------------------------------------------------------------- | -------------------------------------------------------- |
        | [`#!python add_disorder(sublattice [, ...])`][disorder-add_disorder]                  | Add the disorder to the lattice.                         |
        | [`#!python add_local_disorder(sublattice_name [, ...])`][disorder-add_local_disorder] | Internal function to add the disorder to the positions.  |
    
    :   !!! declaration-function "<span id='disorder-add_disorder'>*function* `#!python add_disorder(sublattice, dis_type, mean_value, standard_deviation=0.)`</span>"
            
            
        :   Add the disorder to the lattice.
            
             **Parameters**

             :  | Parameter                                                                  | Description                                                                                                                                                                                              |
                |----------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python sublattice`:*`#!python str` or `#!python list(str)`*             | Name of the sublattice to give a unique value.                                                                                                                                                           |
                | `#!python dis_type`:*`#!python int` or `#!python list(int)`*               | The type of disorder to apply, possible values are `#!python "Gaussian"`, `#!python "Uniform"`, `#!python "Deterministic"`,  `#!python "gaussian"`,  `#!python "uniform"` or `#!python "deterministic"`. |
                | `#!python dis_mean_valuetype`:*`#!python float` or `#!python list(float)`* | Mean value of the deformation.                                                                                                                                                                           |
                | `#!python standard_deviation`:*`#!python float` or `#!python list(float)`* | Standard deviation of the deformation.                                                                                                                                                                   |

    :   !!! declaration-function "<span id='disorder-add_local_disorder'>*function*`#!python add_local_disorder(sublattice_name, dis_type, mean_value, standard_deviation)`</span>"
            
            
        :   Internal function to add the disorder to the positions.
            
             **Parameters**

             :  | Parameter                                              | Description                                                                                                                                                                                              |
                |--------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python sublattice`:*`#!python list(str)`*           | Name of the sublattices to give a unique value.                                                                                                                                                          |
                | `#!python dis_type`:*`#!python list(int)`*             | The type of disorder to apply, possible values are `#!python "Gaussian"`, `#!python "Uniform"`, `#!python "Deterministic"`,  `#!python "gaussian"`,  `#!python "uniform"` or `#!python "deterministic"`. |
                | `#!python dis_mean_valuetype`:*`#!python list(float)`* | Mean value of the deformation.                                                                                                                                                                           |
                | `#!python standard_deviation`:*`#!python list(float)`* | Standard deviation of the deformation.                                                                                                                                                                   |

## Modification
!!! declaration-class "*class* `#!python kite.Modification(magnetic_field=None, flux=None)`"


:   Class that modifies the initially built [`#!python pb.Lattice`][lattice] with a [magnetic field][magnetic-field].

:   **Parameters**
    : <span id="modification-par-magnetic_field">`#!python magnetic_field`: *`#!python float`*</span>
        : Add the magnetic field to the lattice. The field will point along the second primitive lattice vector of the lattice. The magnetic field is in units of $Tesla$, if the [`#!python pb.Lattice`][lattice] is in units of $nm$. The magnetic field is rounded down to the nearest flux quantum.
    : <span id="modification-par-flux">`#!python flux`: *`#!python float`*</span>
        : Add the magnetic flux to the lattice.

:   **Attributes**
    :   | Attribute                                                                                      | Description                                                                                                                                      |
        | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
        | <span id="modification-atr-magnetic_field">`#!python magnetic_field`:*`#!python float`*</span> | The added magnetic field to the lattice.                                                                                                         |
        | <span id="modification-atr-flux">`#!python flux`:*`#!python float`*</span>                     | The added magnetic flux to the lattice. *This is **not** the exact value used in the calculation, but the value added using the parameter above. |

## Configuration
!!! declaration-class "*class* `#!python kite.Configuration(divisions=(1, 1, 1), length=(1, 1, 1), boundaries=('open', 'open', 'open'), is_complex=False, precision=1, spectrum_range=None, angles=(0,0,0), custom_local=False, custom_local_print=False)`"
    
     
:   Define the basic parameters used in the calculation

:   **Parameters**
    : <span id="configuration-divisions">`#!python divisions`: *`#!python int` or `#!python tuple(int, int)` or `#!python tuple(int, int, int)`*</span>
        : Number of decomposition parts of the system.
    : <span id="configuration-length">`#!python length`: *`#!python int` or `#!python tuple(int, int)` or `#!python tuple(int, int, int)`*</span>
        : Number of unit cells in each direction.
    : <span id="configuration-boundaries">`#!python boundaries`: *`#!python str` or `#!python tuple(str, str)` or `#!python tuple(str, str, str)`*</span>
        : Boundary conditions per direction. Possible values are `#!python "periodic"`, `#!python "open"`,
          `#!python "twisted"`*(this option needs the extra argument `#!python angles` and `#!python "random"`.
          See [Settings: Boundaries][settings-boundaries] for what each mode physically means and when
          to choose one over another.
    : <span id="configuration-is_complex">`#!python is_complex`: *`#!python bool`*</span>
        : Boolean that reflects whether the type of Hamiltonian is complex or not.
    : <span id="configuration-precision">`#!python precision`: *`#!python int`*</span>
        : Integer which defines the precision of the number used in the calculation,
          `#!python float` (`#!python 0`), `#!python double` (`#!python 1`), `#!python long double` (`#!python 2`).
    : <span id="configuration-spectrum_range">`#!python spectrum_range`: *`#!python tuple(float, float)`*</span>
        : Energy scale which defines the scaling factor of all the energy related parameters.
          The scaling is done automatically in the background after this definition.
          If the term is not specified, a rough estimate of the bounds is found.
            
            !!! Warning

                Automatic scaling can lead to segmentation-errors due to an error in [pybinding].

    : <span id="configuration-angles">`#!python angles`: *`#!python float` or `#!python tuple(float, float)` or `#!python tuple(float, float, float)`*</span>
        : The angles used for the twisted boundary conditions when `#!python boundary="twist"` is selected.
          The values of `#!python angle`must be in the interval $[0, M \cdot 2 \pi]$.
    : <span id="configuration-custom_local">`#!python custom_local`: *`#!python bool`*</span>
        : Boolean that reflects whether the calculation should use the user-defined local potential.
    : <span id="configuration-custom_local_print">`#!python custom_local_print`: *`#!python bool`*</span>
        : Boolean that reflects whether the calculation should use output the values for the local potential for the various sublattices of the [`#!python pb.Lattice`][lattice].

:   **Attributes**

    :   | Attribute                                                                                                                                                 | Description                                                                                                                                                                                                                                                                                                                                                   |
        |-----------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        | <span id="configuration-energy_scale">`#!python energy_scale`:*`#!python float`*</span>                                                                   | Returns the energy scale of the hopping parameters.                                                                                                                                                                                                                                                                                                           |
        | <span id="configuration-energy_shift">`#!python energy_shift`:*`#!python float`*</span>                                                                   | Returns the energy shift of the hopping parameters around which the spectrum is centered.                                                                                                                                                                                                                                                                     |
        | <span id="configuration-comp">`#!python comp`:*`#!python int`*</span>                                                                                     | Returns `#!python 0` if hamiltonian is real and `#!python 1` elsewise.                                                                                                                                                                                                                                                                                        |
        | <span id="configuration-prec">`#!python prec`:*`#!python int`*</span>                                                                                     | Returns `#!python 0`, `#!python 1`, `#!python 2` if precision if `#!python float`, `#!python double`, and `#!python long double` respectively.                                                                                                                                                                                                                |
        | <span id="configuration-div">`#!python div`:*`#!python int`*</span>                                                                                       | Returns the number of decomposed elements of matrix in $x$, $y$ and/or $z$ direction. Their product gives the total number of threads spawn.                                                                                                                                                                                                                  |
        | <span id="configuration-bound">`#!python bound`:*`#!python tuple(array_like, array_like)`*</span>                                                         | Returns the boundary conditions in each direction, the first argument describes the boundary conditions for the various dimensions with `#!python 0` - open boundary condtions, `#!python 1` - periodic or twisted boundary conditions, `#!python 2` - random boundary conditions, the second gives the angle used if `#!python boundary="twist"` was chosen. |
        | <span id="configuration-leng">`#!python leng`:*`#!python array_like`*</span>                                                                              | Return the number of unit cell repetitions in each direction.                                                                                                                                                                                                                                                                                                 |
        | <span id="configuration-type">`#!python type`:*`#!python np.float32` or `#!python np.float64` or `#!python np.float128` or `#!python np.float256`*</span> | Return the type of the Hamiltonian complex or real, and float, double or long double.                                                                                                                                                                                                                                                                         |
        | <span id="configuration-custom_pot">`#!python custom_pot`:*`#!python bool`*</span>                                                                        | Return custom potential flag.                                                                                                                                                                                                                                                                                                                                 |
        | <span id="configuration-print_custom_pot">`#!python print_custom_pot`:*`#!python bool`*</span>                                                            | Return print custom potential flag.                                                                                                                                                                                                                                                                                                                           |


:   **Methods**
    :   | Method                                         | Description                                                                     |
        | ----------------------------------------------- | ------------------------------------------------------------------------------- |
        | [`#!python set_type()`][configuration-set_type] | Internal function to determine the precision to be used during the calculation. |

    :   !!! declaration-function "<span id="configuration-set_type">*function* `#!python set_type()`</span>"
            
            
        :   Internal function to determine the precision to be used during the calculation.

## Calculation
!!! declaration-class "*class* `#!python kite.Calculation(configuration=None)`"
    
     
:   Class that contains the required target functions to calculate in later steps.

:   **Parameters**
    : <span id="calculation-configuration">`#!python configuration`: *[`#!python kite.Configuration`][configuration]*</span>
        : A [`#!python kite.Configuration`][configuration] that contains the settings for the calculation.

:   **Attributes**
    :   | Attribute                                                                                                                        | Description                                                                                                                 |
        |----------------------------------------------------------------------------------------------------------------------------------| --------------------------------------------------------------------------------------------------------------------------- |
        | <span id="calculation-get_dos">`#!python get_dos`:*`#!python dict`*</span>                                                       | Returns the requested DOS functions.                                                                                        |
        | <span id="calculation-get_ldos">`#!python get_ldos`:*`#!python dict`*</span>                                                     | Returns the requested LDOS functions.                                                                                       |
        | <span id="calculation-get_ldos_map">`#!python get_ldos_map`:*`#!python dict`*</span>                                             | Returns the requested LDOS-map functions.                                                                                   |
        | <span id="calculation-get_spectral_map">`#!python get_spectral_map`:*`#!python dict`*</span>                                     | Returns the requested spectral-map functions.                                                                               |
        | <span id="calculation-get_arpes">`#!python get_arpes`:*`#!python dict`*</span>                                                   | Returns the requested ARPES functions.                                                                                      |
        | <span id="calculation-get_gaussian_wave_packet">`#!python get_gaussian_wave_packet`:*`#!python dict`*</span>                     | Returns the requested wave packet time evolution function, with a gaussian wavepacket mutiplied with different plane waves. |
        | <span id="calculation-get_localized_wave_packet">`#!python get_localized_wave_packet`:*`#!python dict`*</span>                   | Returns the requested localized wave packet time evolution function, with spectral filtering.                              |
        | <span id="calculation-get_conductivity_dc">`#!python get_conductivity_dc`:*`#!python dict`*</span>                               | Returns the requested DC conductivity functions.                                                                            |
        | <span id="calculation-get_conductivity_optical">`#!python get_conductivity_optical`:*`#!python dict`*</span>                     | Returns the requested optical conductivity functions.                                                                       |
        | <span id="calculation-get_conductivity_optical_nonlinear">`#!python get_conductivity_optical_nonlinear`:*`#!python dict`*</span> | Returns the requested nonlinear optical conductivity functions.                                                             |
        | <span id="calculation-get_singleshot_conductivity_dc">`#!python get_singleshot_conductivity_dc`:*`#!python dict`*</span>         | Returns the requested singleshot DC conductivity functions.                                                                 |
        | <span id="calculation-get_custom_one">`#!python get_custom_one`:*`#!python dict`*</span>                                         | Returns the requested rank-one custom-operator trace functions.                                                             |
        | <span id="calculation-get_custom_one_local">`#!python get_custom_one_local`:*`#!python dict`*</span>                             | Returns the requested local (position/energy-resolved) rank-one custom-operator functions.                                  |
        | <span id="calculation-get_custom_two">`#!python get_custom_two`:*`#!python dict`*</span>                                         | Returns the requested rank-two custom-operator trace functions.                                                             |
        | <span id="calculation-get_custom_two_local">`#!python get_custom_two_local`:*`#!python dict`*</span>                             | Returns the requested local rank-two custom-operator functions.                                                             |
        | <span id="calculation-get_custom_ss_two">`#!python get_custom_ss_two`:*`#!python dict`*</span>                                   | Returns the requested single-shot rank-two custom-operator functions.                                                       |


:   **Methods**
    :   | Method                                                                                         | Description                                                                 |
        |------------------------------------------------------------------------------------------------| --------------------------------------------------------------------------- |
        | [`#!python dos(num_points, num_moments, [, ...]`)][calculation-dos]                            | Calculate the density of states as a function of energy.                    |
        | [`#!python ldos(energy, num_moments, [, ...])`][calculation-ldos]                              | Calculate the local density of states as a function of energy.              |
        | [`#!python ldos_map(energy_, sigma_, vectors_ [, ...])`][calculation-ldos_map]                  | Calculate a real-space map of the local density of states at one energy.    |
        | [`#!python spectral_map(energy_, sigma_, vectors_ [, ...])`][calculation-spectral_map]          | Calculate a momentum-space (k-resolved) map of the spectral function at one energy. |
        | [`#!python arpes(k_vector, weight [, ...])`][calculation-arpes]                                | Calculate the spectral contribution for given k-points and weights.         |
        | [`#!python gaussian_wave_packet(num_points [, ...])`][calculation-gaussian_wave_packet]        | Calculate the time evolution function of a wave packet.                     |
        | [`#!python localized_wave_packet(time, num_measures, initial_pos [, ...])`][calculation-localized_wave_packet] | Calculate the time evolution of a localized/gaussian wave packet with spectral filtering. |
        | [`#!python conductivity_dc(direction, [, ...])`][calculation-conductivity_dc]                  | Calculate the DC conductivity for a given direction.                        |
        | [`#!python conductivity_optical(direction, [, ...])`][calculation-conductivity_optical]        | Calculate optical conductivity for a given direction.                       |
        | [`#!python conductivity_optical_nonlinear([...])`][calculation-conductivity_optical_nonlinear] | Calculate nonlinear optical conductivity for a given direction.             |
        | [`#!python singleshot_conductivity_dc(energy, [...])`][calculation-singleshot_conductivity_dc] | Calculate the DC conductivity using KITEx for a given direction and energy. |
        | [`#!python add_orbital_index(label_, idx_)`][calculation-add_orbital_index]                    | Register a name for an orbital index, for use in custom operator strings.   |
        | [`#!python add_orbital_coupling(start_, last_, c_, label_)`][calculation-add_orbital_coupling] | Define one matrix element of a custom orbital-space operator.               |
        | [`#!python custom_one(stream_, num_random_, num_disorder_)`][calculation-custom_one]           | Calculate the rank-one custom-operator trace `#!python Tr[Tn(H)·J]`.        |
        | [`#!python custom_one_local(stream_, energy_, position_, sublattice_ [, ...])`][calculation-custom_one_local] | Calculate the rank-one custom-operator trace at a chosen energy and position. |
        | [`#!python custom_two(stream_, num_random_, num_disorder_, num_points_, temperature_)`][calculation-custom_two] | Calculate the rank-two custom-operator trace `#!python Tr[Tn(H)·A·Tm(H)·B]`. |
        | [`#!python custom_two_local(stream_, positions_)`][calculation-custom_two_local]               | Calculate the rank-two custom-operator trace at chosen positions.           |
        | [`#!python custom_singleshot_two(stream_, num_random_, num_disorder_ [, ...])`][calculation-custom_singleshot_two] | Calculate the rank-two custom-operator trace at chosen energies, using KITEx (single-shot method). |

    :   !!! declaration-function "<span id="calculation-dos">*function* `#!python dos(num_points, num_moments, num_random, num_disorder=1)`</span>"
            
            
        :   Calculate the density of states as a function of energy.
            
            **Parameters**

            :   | Parameter                                | Description                                                                     |
                |------------------------------------------|---------------------------------------------------------------------------------|
                | `#!python num_points`:*`#!python int`*   | Number of energy point inside the spectrum at which the DOS will be calculated. |
                | `#!python num_moments`:*`#!python int`*  | Number of polynomials in the Chebyshev expansion.                               |
                | `#!python num_random`:*`#!python int`*   | Number of random vectors to use for the stochastic evaluation of trace.         |
                | `#!python num_disorder`:*`#!python int`* | Number of different disorder realisations.                                      |

    
    :   !!! declaration-function "<span id="calculation-ldos">*function*`#!python ldos(energy, num_moments, position, sublattice, num_disorder=1)`</span>"
            
            
        :   Calculate the local density of states as a function of energy. 
            
             **Parameters**

            :   | Parameter                                                   | Description                                                        |
                |-------------------------------------------------------------|--------------------------------------------------------------------|
                | `#!python energy`:*`#!python array_like`*                   | List of energy points at which the LDOS will be calculated.        |
                | `#!python num_moments`:*`#!python int`*                     | Number of polynomials in the Chebyshev expansion.                  |
                | `#!python position`:*`#!python int`*                        | Relative index of the unit cell where the LDOS will be calculated. |
                | `#!python sublattice`:*`#!python list`*                     | Name of the sublattice at which the LDOS will be calculated.       |
                | `#!python num_disorder`:*`#!python str` or `#!python list`* | Number of different disorder realisations.                         |

    :   !!! declaration-function "<span id="calculation-ldos_map">*function*`#!python ldos_map(energy_, sigma_, vectors_, coef="gaussian")`</span>"


        :   Calculate a full real-space map of the local density of states at one target energy.

            KITE evaluates this at every lattice site simultaneously with a stochastic (random-phase-vector)
            Chebyshev/KPM estimator combined with a Gaussian- or window-broadened energy filter $f(H)$.
            Up to the overall normalization, the stored quantity is $\mathbb{E}[\text{map}_r] = \text{factor}\cdot\langle r|f(H)^2|r\rangle$
            — a filtered-random-state local-intensity estimator, *not* simply $\langle r|f(H)|r\rangle$. This is a standard
            stochastic-trace technique; see the review by Weiße, Wellein, Alvermann, and Fehske.[^1]

            **Parameters**

            :   | Parameter                                   | Description                                                                                          |
                |-----------------------------------------------|-------------------------------------------------------------------------------------------------------|
                | `#!python energy_`:*`#!python float`*         | Target energy where the map is evaluated.                                                             |
                | `#!python sigma_`:*`#!python float`*          | Width of the Gaussian (or window) that approximates the Dirac delta.                                  |
                | `#!python vectors_`:*`#!python int`*          | Number of independent random-phase vectors to average over (analogous to `#!python num_random` for [`#!python dos()`][calculation-dos] — *not* a number of k-points or spatial points). |
                | `#!python coef`:*`#!python str`*              | Choice of energy-filter kernel, either `#!python "gaussian"` or `#!python "window"`.                  |

            !!! Note "Disorder is not ensemble-averaged"

                Unlike [`#!python dos()`][calculation-dos] and [`#!python arpes()`][calculation-arpes], `#!python ldos_map()`
                has no `#!python num_disorder` parameter: a single disorder realisation is generated per call. This is a
                deliberate design choice appropriate for visualizing the real-space structure of *one* disorder
                realization (e.g. impurity resonances), not an oversight — but the output of a single call should not be
                interpreted as a disorder-ensemble average.

            !!! Warning "`is_complex=True` is mandatory, even for real-valued hoppings"

                `#!python ldos_map()` (and `#!python spectral_map()`) are only compiled into KITEx for the complex
                Hamiltonian instantiation (`#!cpp if constexpr (is_tt<std::complex, T>::value)` in
                `Src/Simulation/SimulationLDoS.cpp`). With `#!python is_complex=False`, that branch is compiled away
                *entirely* — KITEx runs, prints "Calculating LDoS. Done.", and exits normally, but **no**
                `#!python /Calculation/ldos_map/Map` dataset is ever written. This is a silent no-op, not an error:
                verified directly by running with `#!python is_complex=False` and finding the expected dataset simply
                absent from the output file. Always set `#!python is_complex=True` for these two methods, even if
                every hopping in your lattice is real.

            !!! Info "Example"

                `#!python examples/piflux_ldos_map.py` maps the LDOS at the Dirac-point energy ($E=0$) around a
                single vacancy in a $\pi$-flux square lattice — the LDOS vanishes exactly at the vacancy site and
                shows a decaying, oscillatory enhancement around it. See the [in-depth write-up][markov_maps_example]
                for the full explanation and plots. This is the reference method for stochastic real-space/momentum-space
                maps with a rigorous, Markov-inequality-bounded per-site sampling error.[^3]

    :   !!! declaration-function "<span id="calculation-spectral_map">*function*`#!python spectral_map(energy_, sigma_, vectors_, coef="gaussian")`</span>"


        :   Calculate a full momentum-space (k-resolved) map of the spectral function at one target energy — the
            k-resolved analogue of [`#!python ldos_map()`][calculation-ldos_map]. It uses the same random-phase/Chebyshev-filter
            machinery, but wraps the filtering loop in a unitary FFT that transforms between the real-space (site) and
            momentum-space (k) bases — a proper Bloch transform that includes intra-cell sublattice phases, not a naive
            plane DFT — before and after applying $H$: $\mathbb{E}[\text{map}_k] = \text{factor}\cdot\langle k|f(H)^2|k\rangle$.

            **Parameters**

            :   | Parameter                                   | Description                                                                                          |
                |-----------------------------------------------|-------------------------------------------------------------------------------------------------------|
                | `#!python energy_`:*`#!python float`*         | Target energy where the map is evaluated.                                                             |
                | `#!python sigma_`:*`#!python float`*          | Width of the Gaussian (or window) that approximates the Dirac delta.                                  |
                | `#!python vectors_`:*`#!python int`*          | Number of independent random-phase vectors to average over.                                           |
                | `#!python coef`:*`#!python str`*              | Choice of energy-filter kernel, either `#!python "gaussian"` or `#!python "window"`.                  |

            !!! Info "Relation to `arpes()`"

                `#!python spectral_map()` is conceptually the "all-k-at-once" stochastic sibling of
                [`#!python arpes()`][calculation-arpes], which instead evaluates exact Chebyshev moments at a handful of
                user-chosen k-points via a genuine plane-wave state (no FFT, no randomness). The FFT in `spectral_map()`
                only does non-trivial work because real-space disorder breaks translation symmetry inside the sandwiched
                Chebyshev recursion — for a clean/periodic system it would trivially return the same map at every k-point
                regardless of averaging.

            !!! Warning "Disorder and boundary conditions are not ensemble-averaged consistently"

                As with [`#!python ldos_map()`][calculation-ldos_map], disorder is frozen per call (no `#!python num_disorder`).
                In addition, if you configure random or twisted boundary conditions (see [`#!python kite.Configuration`][configuration]),
                each realization samples a different, randomly-twisted boundary condition while the FFT's k-grid is fixed
                from the untwisted lattice size — averaging realizations under that combination is not accounted for by
                the normalization and may not do what you expect. Fixed (non-twisted) boundaries are recommended when using
                `#!python spectral_map()` unless you specifically understand this interaction.

            !!! Warning "The FFT runs over every axis, regardless of its declared boundary condition"

                `Src/FFT/FFT.cpp` performs the momentum-space transform over *all* axes unconditionally, even an
                axis declared `#!python "open"` (a real, hard-wall boundary where momentum isn't actually conserved).
                Verified directly in `#!python examples/weyl_spectral_map.py`: with `#!python boundaries=["open", "random", "random"]`,
                the spectral weight along the two periodic/twisted axes is sharply peaked, while along the open axis
                it is smeared and nearly flat — a real, physically-expected consequence of transforming an axis with
                no translational symmetry, not a bug, but easy to misread as a broken/noisy result if you don't expect it.

            !!! Info "Example"

                `#!python examples/weyl_spectral_map.py` computes $A(\mathbf{k}, E{=}0.7)$ for a 3D Weyl semimetal and
                shows two clean ring contours — constant-energy cuts through the two Weyl nodes' linearly-dispersing
                cones. See the [in-depth write-up][markov_maps_example] for the full explanation and plots.

            See the base KITE method paper[^2] for the general LDOS/ARPES-of-disordered-materials methodology
            demonstrated by this family of target functions; the specific real-space↔k-space basis-change construction
            used internally by `#!python spectral_map()` does not have a separately located dedicated reference.

    :   !!! declaration-function "<span id="calculation-arpes">*function*`#!python arpes(k_vector, weight, num_moments, num_disorder=1)`</span>"
            
            
        :   Calculate the spectral contribution for given k-points and weights.
            
            **Parameters**

            :   | Parameter                                                   | Description                                                                                                   |
                |-------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
                | `#!python k_vector`:*`#!python array_like`*                 | List of K points with respect to reciprocal vectors b0 and b1 at which the band structure will be calculated. |
                | `#!python weight`:*`#!python array_like`*                   | List of orbital weights used for ARPES.                                                                       |
                | `#!python num_moments`:*`#!python int`*                     | Number of polynomials in the Chebyshev expansion.                                                             |
                | `#!python num_disorder`:*`#!python int`*                    | Number of different disorder realisations.                                                                    |

    
    :   !!! declaration-function "<span id="calculation-gaussian_wave_packet">*function*`#!python gaussian_wave_packet(num_points, num_moments, timestep, k_vector, spinor, width, mean_value, num_disorder=1, operators=None, probing_point=0)`</span>"
            
            
        :   Calculate the time evolution function of a wave packet.
            
            **Parameters**

            :   | Parameter                                                          | Description                                                                                                     |
                |--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
                | `#!python num_points`:*`#!python int`*                             | Number of time points for the time evolution.                                                                   |
                | `#!python num_moments`:*`#!python int`*                            | Number of polynomials in the Chebyshev expansion.                                                               |
                | `#!python timestep`:*`#!python float`*                             | Timestep for calculation of time evolution.                                                                     |
                | `#!python k_vector`:*`#!python array_like`*                        | Different wave vectors, components corresponding to vectors b0 and b1.                                          |
                | `#!python spinor`:*`#!python array_like`*                          | Spinors for each of the k vectors.                                                                              |
                | `#!python width`:*`#!python float`*                                | Width of the gaussian.                                                                                          |
                | `#!python mean_value`:*`#!python tuple(float, float)`*             | Mean value of the gaussian envelope.                                                                            |
                | `#!python num_disorder`:*`#!python int`*                           | Number of different disorder realisations.                                                                      |
                | `#!python operators`:*`#!python list[str]`*                       | Labels of operators registered via [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling], whose expectation value is tracked at every timestep (e.g. `#!python ['l0', 'l1', 'l2']` for a 3-component spin, orbital-angular-momentum, or quadrupole operator set). |
                | `#!python probing_point`:*`#!python int` or `#!python array_like`* | Forward probing point, defined with x, y coordinate were the wavepacket will be checked at different timesteps. |

            **Tracked operators**

            :   `#!python gaussian_wave_packet()` has no built-in notion of spin: any operator, including spin, is
                tracked through the same [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling]
                mechanism used by [`#!python custom_one()`][calculation-custom_one]/[`#!python custom_two()`][calculation-custom_two]
                — including its **on-site-only restriction** (see the warning under
                [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling]): a tracked operator must be
                expressible as a single matrix acting on one site's orbitals, identical at every site. Spin, orbital
                angular momentum, and orbital quadrupole operators satisfy this; an operator that couples different
                sites cannot be tracked this way. Each requested label's expectation value
                $\langle\psi(t)|\hat O|\psi(t)\rangle$ is written to `#!python /Calculation/gaussian_wave_packet/<label>`
                in the output file. Unlike `#!python custom_one()`'s `#!python stream_` tokens, these labels are read
                back by an explicit ordered list, not the
                single-digit `#!python "l0"`–`#!python "l9"` scheme described in the warning above — so this path
                does not share that 10-operator limit.

    :   !!! declaration-function "<span id="calculation-localized_wave_packet">*function*`#!python localized_wave_packet(time, num_measures, initial_pos, num_moments=0, width=-1., energy_window=[0.,0.], initial_wavevector=None, probes=None, sample_start=-1, sample_L=-1)`</span>"


        :   Calculate the time evolution of a localized (point-like) or Gaussian-envelope wave packet, with optional
            spectral filtering. This is distinct from [`#!python gaussian_wave_packet()`][calculation-gaussian_wave_packet]:
            the initial state is dispatched between a point-like or a Gaussian-envelope wave packet depending on
            whether `#!python width > 0`, it can be band/energy-filtered to `#!python energy_window` before propagation,
            and it reports richer observables (spatial spread statistics, return probability, propagator amplitude,
            spectral moments, and transmission weights at chosen probe positions) than the plain spin expectation
            values returned by `#!python gaussian_wave_packet()`.

            **Parameters**

            :   | Parameter                                              | Description                                                                                                                          |
                |----------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python time`:*`#!python float`*                       | Total dimensionless time over which the wave packet is evolved. Together with `#!python num_measures`, this sets the time step between measurements. |
                | `#!python num_measures`:*`#!python int`*                 | Number of times the evolved state is probed.                                                                                          |
                | `#!python initial_pos`:*`#!python array_like`*           | Initial position of the localized wave packet in lattice coordinates, `#!python [n_x, n_y[, n_z], n_orbital]`.                       |
                | `#!python num_moments`:*`#!python int`*                  | Number of Chebyshev moments to calculate for the spectral function of the propagated state.                                          |
                | `#!python width`:*`#!python float`*                      | Width of the Gaussian envelope in real space. If `#!python width <= 0` a point-like (not Gaussian-enveloped) wave packet is used.     |
                | `#!python energy_window`:*`#!python array_like`*         | The localized packet is filtered to keep only eigenstates with energy inside this window.                                            |
                | `#!python initial_wavevector`:*`#!python array_like`*    | Wavevector the Gaussian envelope is centered around, in reciprocal lattice coordinates.                                               |
                | `#!python probes`:*`#!python array_like`*                | List of positions, in lattice coordinates, where the propagator/transmission weights are evaluated.                                  |
                | `#!python sample_start`:*`#!python int`*                 | Start (in the x-direction) of a disordered sample region embedded between clean leads, for transport/localization-length studies.     |
                | `#!python sample_L`:*`#!python int`*                     | Length of the disordered sample region described by `#!python sample_start`.                                                          |

            !!! Note "No worked example yet"

                There is currently no example script in `#!bash examples/` that uses `#!python localized_wave_packet()`.
                See the C++ implementation (`Src/Simulation/SimulationLocalizedWavePacket.cpp`) and the parameter reference
                above for now.

    :   !!! declaration-function "<span id="calculation-conductivity_dc">*function*`#!python conductivity_dc(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0)`</span>"
            
            
        :   Calculate the DC conductivity for a given direction.

            !!! Info "Output units: e²/h"

                The reconstructed conductivity (via [`#!bash KITE-tools --CondDC`][kitetools])
                is in units of $e^2/h$ directly — **not** $e^2/\hbar$, $e^2/(\pi h)$, or any other
                common convention. Verified empirically (not just asserted): the Hall conductivity
                $\sigma_{xy}$ of a Haldane Chern insulator ($C=1$,
                `#!python examples/dos_dccond_haldane.py`'s own parameters, including its Anderson
                disorder) reconstructs to $\sigma_{xy}\approx1.01$–$1.02$ at mid-gap — a topological
                plateau quantized at exactly $C\cdot e^2/h=1\,e^2/h$, insensitive to broadening,
                disorder, or moment-count details, unlike a semiclassical Drude-regime check would
                be. The same convention applies to [`#!python conductivity_optical()`][calculation-conductivity_optical]
                and the [`#!python custom_two()`][calculation-custom_two] family (all built from the
                same underlying Kubo-Bastin/Gamma2D machinery).

            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                      |
                |-------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"`, `#!python "zz"`, `#!python "xy"`, `#!python "xz"`, `#!python "yx"`, `#!python "yz"`, `#!python "zx"`, `#!python "zy"`.     |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                  |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                          |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                       |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy. |

    
    
    
    :   !!! declaration-function "<span id="calculation-conductivity_optical">*function*`#!python conductivity_optical(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0)`</span>"
            
            
        :   Calculate optical conductivity for a given direction. Same output-units convention as
            [`#!python conductivity_dc()`][calculation-conductivity_dc] ($e^2/h$).
                        
            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                      |
                |-------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"`, `#!python "zz"`, `#!python "xy"`, `#!python "xz"`, `#!python "yx"`, `#!python "yz"`, `#!python "zx"`, `#!python "zy"`.     |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                  |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                          |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                       |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy. |

    
    :   !!! declaration-function "<span id="calculation-conductivity_optical_nonlinear">*function*`#!python conductivity_optical_nonlinear(direction, num_points, num_moments, num_random, num_disorder=1, temperature=0, special=0)`</span>"
            
            
        :   Calculate nonlinear optical conductivity for a given direction.
            
            **Parameters**

            :   | Parameter                                 | Description                                                                                                                                                                                                                                                                  |
                |-------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
                | `#!python direction`:*`#!python str`*     | Direction in $xyz$-coordinates along which the conductivity is calculated, supports all the combinations of the directions `#!python "x"`, `#!python "y"` and `#!python "z"` with length 3 like `#!python "xxx"`, `#!python "zzz"`, `#!python "xxy"`, `#!python "xxz"` etc.  |
                | `#!python num_points`:*`#!python int`*    | Number of energy point inside the spectrum at which the DOS will be calculated.                                                                                                                                                                                              |
                | `#!python num_moments`:*`#!python int`*   | Number of polynomials in the Chebyshev expansion.                                                                                                                                                                                                                            |
                | `#!python num_random`:*`#!python int`*    | Number of random vectors to use for the stochastic evaluation of trace.                                                                                                                                                                                                      |
                | `#!python num_disorder`:*`#!python int`*  | Number of different disorder realisations.                                                                                                                                                                                                                                   |
                | `#!python temperature`:*`#!python float`* | Value of the temperature at which we calculate the response. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitraty units, specify the quantity $K_B \cdot T$, which has units of energy.                             |
                | `#!python special`:*`#!python int`*       | Optional, a parameter that can simplify the calculation for some materials.                                                                                                                                                                                                  |
    
    
    :   !!! declaration-function "<span id="calculation-singleshot_conductivity_dc">*function*`#!python singleshot_conductivity_dc(energy, direction, eta, num_moments, num_random, num_disorder=1, preserve_disorder=False)`</span>"
            
            
        :   Calculate the DC conductivity using KITEx for a given direction and energy.
            
            !!! Info "Processing the output of `#!python singleshot_conductivity_dc()`"

                `#!python singleshot_conductivity_dc()` works different from the other target-functions in that a
                single run with [KITEx][kitex] is sufficient. The results don't have to be processed by
                [KITE-tools][kitetools].
                As such, the results are already available in the [HDF5]-file. You can extract the results from the
                [HDF5]-file [as explained in the tutorial][tutorial-hdf5].

                There is also [a script][kitetools-output] that automates this process with, `#!python "output.h5"`
                the name of the [HDF5]-file processed by [KITEx][kitex]:
                
                ``` bash
                python3 process_single_shot.py output.h5
                ``` 
                
                this will result in a data-file `#!python "output.dat"`, as explained in the
                [API of KITE-tools][kitetools-output].                

            **Parameters**

            :   | Parameter                                                     | Description                                                                                                                               |
                |---------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python energy`:*`#!python array_like` or `#!python float`* | Array or a single value of energies at which `#!python singleshot_conductivity_dc` will be calculated.                                    |
                | `#!python direction`:*`#!python str`*                         | Direction in $xyz$-coordinates along which the conductivity is calculated, supports `#!python "xx"`, `#!python "yy"` and `#!python "zz"`. |
                | `#!python eta`:*`#!python int`*                               | Parameter that affects the broadening of the kernel function.                                                                             |
                | `#!python num_moments`:*`#!python int`*                       | Number of polynomials in the Chebyshev expansion.                                                                                         |
                | `#!python num_random`:*`#!python int`*                        | Number of random vectors to use for the stochastic evaluation of trace.                                                                   |
                | `#!python num_disorder`:*`#!python int`*                      | Number of different disorder realisations.                                                                                                |
                | `#!python preserve_disorder`:*`#!python bool`*                | Optional.                                                                                                                                 |

    :   !!! declaration-function "<span id="calculation-add_orbital_index">*function*`#!python add_orbital_index(label_, idx_)`</span>"


        :   Register a name for an orbital index, for later reference inside a custom operator string
            (see [`#!python custom_one()`][calculation-custom_one] below). Must be called once per orbital before
            using [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling].

            **Parameters**

            :   | Parameter                          | Description                                                                                          |
                |--------------------------------------|--------------------------------------------------------------------------------------------------------|
                | `#!python label_`:*`#!python str`*   | Name used to refer to this orbital when calling [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling]. |
                | `#!python idx_`:*`#!python int`*     | Index (0-based) assigned to this orbital inside the custom orbital-space operator matrices.            |

    :   !!! declaration-function "<span id="calculation-add_orbital_coupling">*function*`#!python add_orbital_coupling(start_, last_, c_, label_)`</span>"


        :   Set one matrix element of a named, custom orbital-space operator matrix $M_{\texttt{label\_}}$, for use as an
            `#!python "l"`-type token inside a custom operator string (see below). Internally this sets
            $M_{\texttt{label\_}}[\text{idx}(\texttt{last\_}),\,\text{idx}(\texttt{start\_})] = \texttt{c\_}$, where
            `#!python idx(...)` is the mapping registered with [`#!python add_orbital_index()`][calculation-add_orbital_index].
            Calling this repeatedly with the same `#!python label_` accumulates further matrix elements into the same
            operator matrix (which is otherwise initialized to zero).

            **Parameters**

            :   | Parameter                            | Description                                                                                                        |
                |----------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
                | `#!python start_`:*`#!python str`*     | Name (registered via [`#!python add_orbital_index()`][calculation-add_orbital_index]) of the orbital coupled *from*. |
                | `#!python last_`:*`#!python str`*      | Name of the orbital coupled *to*.                                                                                     |
                | `#!python c_`:*`#!python complex`*     | Value of the matrix element.                                                                                          |
                | `#!python label_`:*`#!python str`*     | Name of this custom orbital operator; must begin with `#!python "l"` (raises `#!python ValueError` otherwise).       |

            !!! Warning "`#!python \"l\"`-type operators are on-site only"

                $M_{\texttt{label\_}}$ is a single orbital-space matrix, applied identically and independently at
                *every* lattice site — it can only mix orbitals belonging to the same site, never couple different
                sites or unit cells. This is why it's a separate mechanism from the `#!python "v"` (velocity) token:
                velocity is intrinsically an inter-site (hopping) operator, computed a completely different way
                (`#!python h.build_velocity()`), not expressible as a per-site matrix. Spin, orbital angular
                momentum, and orbital quadrupole operators are all valid `#!python "l"`-type operators because they
                are genuinely on-site quantities (built entirely from a single site's orbital manifold, e.g. the
                $p_x,p_y,p_z$ orbitals of one atom) — but this mechanism cannot represent an operator that is
                diagonal in *orbital* index yet varies from site to site, nor one with any inter-site matrix element.

            !!! Info "The `kite.custom.Vertex` operator language"
                <span id="calculation-vertex-grammar"></span>

                A [`#!python kite.custom.Vertex(num_moments, stream)`][calculation-custom_one] object (import with
                `#!python from kite import custom`) defines a custom operator $J = \sum_i c_i O_i$ as a weighted sum of
                operator-string products, to be traced against Chebyshev polynomials of $H$. `#!python num_moments` is
                the number of Chebyshev moments computed for this vertex — note that in the shipped example scripts this
                argument is passed as a variable named `#!python pol_A`/`#!python pol_B`, which is *not* a polarization,
                it is a moment count. `#!python stream` is a list of `#!python [coefficient, operator_string]` pairs;
                operators can also be appended afterwards with `#!python vertex.add_operator(coef, operator_string)`.

                Each `#!python operator_string` is a dot-separated chain of tokens, e.g. `#!python "vy.rx"`. The first
                character of each token selects the operator type:

                - `#!python "v" + component` (e.g. `#!python "vy"`) — the velocity operator along that Cartesian direction
                  (`#!python "x"`, `#!python "y"`, or `#!python "z"`).
                - `#!python "r" + component` (e.g. `#!python "rx"`) — the real-space position operator along that direction.
                - `#!python "l" + digit` (e.g. `#!python "l0"`) — a custom orbital-space matrix registered with
                  [`#!python add_orbital_index()`][calculation-add_orbital_index] /
                  [`#!python add_orbital_coupling()`][calculation-add_orbital_coupling].

                A dotted chain composes the tokens as an operator product applied to the state right-to-left as written,
                i.e. `#!python "vy.rx"` applies `#!python "rx"` to the state first and then `#!python "vy"` — equivalent
                to the ordinary matrix product $v_y r_x$ acting on a ket.

                !!! Warning "Only single-digit orbital-operator indices are supported"

                    Because the C++ reader consumes exactly one character after `#!python "l"`, only labels
                    `#!python "l0"`–`#!python "l9"` can be addressed inside an operator string — a hard limit of 10
                    distinct custom orbital operators per model. The token's digit selects the orbital operator
                    *positionally* from the operators exported to the HDF5 file; this was not independently verified
                    end-to-end against the HDF5 read-back order beyond the `#!python "l0", "l1", ...` naming convention,
                    so it is safest to register your custom orbital operators as `#!python "l0"`, `#!python "l1"`,
                    `#!python "l2"`, … in the same order you intend to reference them.

    :   !!! declaration-function "<span id="calculation-custom_one">*function*`#!python custom_one(stream_, num_random_, num_disorder_)`</span>"


        :   Calculate the generalized rank-one KPM trace `#!python Tr[Tn(H)·J]` for an arbitrary custom operator $J$
            built from a [`#!python kite.custom.Vertex`][calculation-custom_one] (see the operator-language box above).

            **Parameters**

            :   | Parameter                                | Description                                                                     |
                |---------------------------------------------|------------------------------------------------------------------------------------|
                | `#!python stream_`:*`#!python kite.custom.Vertex`* | Vertex object defining the operator $J$ as a weighted sum of operator strings.       |
                | `#!python num_random_`:*`#!python int`*     | Number of random vectors to use for the stochastic evaluation of the trace.        |
                | `#!python num_disorder_`:*`#!python int`*   | Number of different disorder realisations.                                         |

            !!! Info "The `#!python \"v\"` (velocity) building block and NumVelocities"

                KITE's raw `#!python "vx"`/`#!python "vy"` DSL token is $[\hat H,\hat r]$ — missing the
                $i/\hbar$ factor of the textbook Hermitian velocity operator. `#!python custom_one()` counts
                how many `#!python "v"`-type tokens appear in `#!python stream_` itself (every term must
                carry the same count, or a `#!python ValueError` is raised) and writes it as `NumVelocities`
                in the output file. This determines whether `#!python Tr[Tn(H)·J]` is physically real (even
                count) or purely imaginary (odd count) — both are genuine signal, not noise — and lets
                [`#!bash KITE-tools --CustomOne`][kite-tools-customone] reconstruct the correct component
                automatically. You do not need to insert a compensating `#!python i` into the vertex
                yourself; see `#!python examples/haldane_orbital_magnetization.py` for a worked example with
                an odd (single) velocity operator.

            !!! Example "Worked example"

                `#!python examples/haldane_orbital_magnetization.py` computes the orbital magnetization of a
                Haldane-model Chern insulator via `#!python custom_one()`, post-processed with
                [`#!bash KITE-tools --CustomOne`][kite-tools-customone].

    :   !!! declaration-function "<span id="calculation-custom_one_local">*function*`#!python custom_one_local(stream_, energy_, position_, sublattice_, num_disorder_=1)`</span>"


        :   The position/energy-resolved (LDOS-like) sibling of [`#!python custom_one()`][calculation-custom_one]: evaluates
            the same rank-one trace at one chosen real-space position and sublattice, as a function of energy, instead of
            as a single lattice-wide trace.

            **Parameters**

            :   | Parameter                                              | Description                                                                     |
                |-----------------------------------------------------------|--------------------------------------------------------------------------------|
                | `#!python stream_`:*`#!python kite.custom.Vertex`*        | Vertex object defining the operator $J$ (see [`#!python custom_one()`][calculation-custom_one]). |
                | `#!python energy_`:*`#!python array_like`*                | List of energy points at which the trace will be calculated.                    |
                | `#!python position_`:*`#!python array_like`*              | Relative index of the unit cell where the trace will be calculated.             |
                | `#!python sublattice_`:*`#!python str` or `#!python list`*| Name of the sublattice at which the trace will be calculated.                  |
                | `#!python num_disorder_`:*`#!python int`*                 | Number of different disorder realisations.                                     |

    :   !!! declaration-function "<span id="calculation-custom_two">*function*`#!python custom_two(stream_, num_random_, num_disorder_, num_points_, temperature_)`</span>"


        :   Calculate the generalized rank-two KPM trace `#!python Tr[Tn(H)·A·Tm(H)·B]` for two independently-defined
            [`#!python kite.custom.Vertex`][calculation-custom_one] operator streams `#!python A`, `#!python B`, passed
            as `#!python stream_=[A, B]`. When `#!python A`, `#!python B` are velocity operators this reduces to
            [`#!python conductivity_dc()`][calculation-conductivity_dc]'s own machinery, so the reconstructed trace
            carries the same $e^2/h$ convention **only when the operators involved have the same physical dimensions
            as velocity** (e.g. the spin-current vertex on the
            [spin Hall page][custom-vertex-example]) — a **bare density operator** vertex (e.g. the
            [Rashba-Edelstein example][rashba-edelstein-example]'s bare spin density) is a genuinely different
            physical observable and needs its own dimensional analysis from the underlying Kubo-Bastin
            formula, not an assumed carry-over of this constant.

            **Parameters**

            :   | Parameter                                         | Description                                                                                                                                                                                                                                      |
                |-------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
                | `#!python stream_`:*`#!python list(kite.custom.Vertex, kite.custom.Vertex)`* | The two vertices `#!python [A, B]` defining the operators $A$ and $B$.                                                                                                                                     |
                | `#!python num_random_`:*`#!python int`*               | Number of random vectors to use for the stochastic evaluation of the trace.                                                                                                                                                                       |
                | `#!python num_disorder_`:*`#!python int`*             | Number of different disorder realisations.                                                                                                                                                                                                        |
                | `#!python num_points_`:*`#!python int`*               | Number of energy points used by post-processing to reconstruct the energy dependence.                                                                                                                                                            |
                | `#!python temperature_`:*`#!python float`*            | Value of the temperature at which the response is calculated. If $eV$ is used as unit for energy, then $k_B\cdot T$ is also in $eV$. To define the temperature in arbitrary units, specify the quantity $K_B \cdot T$, which has units of energy. |

            !!! Info "Example"

                See `#!python examples/kane_mele_spin_hall.py` (and its disorder extension
                `#!python kane_mele_spin_hall_disorder.py`) for a worked example of `#!python custom_two()`,
                walked through in depth in [Custom Vertex Operators][custom-vertex-example].

    :   !!! declaration-function "<span id="calculation-custom_two_local">*function*`#!python custom_two_local(stream_, positions_)`</span>"


        :   The position-resolved sibling of [`#!python custom_two()`][calculation-custom_two]: evaluates the same
            rank-two trace at a list of chosen positions, with no random-vector/disorder averaging, energy-point count,
            or temperature argument.

            **Parameters**

            :   | Parameter                                    | Description                                                                                      |
                |---------------------------------------------------|----------------------------------------------------------------------------------------------------|
                | `#!python stream_`:*`#!python list(kite.custom.Vertex, kite.custom.Vertex)`* | The two vertices `#!python [A, B]` defining the operators $A$ and $B$. |
                | `#!python positions_`:*`#!python array_like`*     | Array of shape `#!python (number_of_positions, D + 1)` listing the positions at which the trace is evaluated ($D$ lattice directions plus a sublattice/orbital index). |

            !!! Note "No worked example yet"

                There is currently no example script that uses `#!python custom_two_local()`.

    :   !!! declaration-function "<span id="calculation-custom_singleshot_two">*function*`#!python custom_singleshot_two(stream_, num_random_, num_disorder_, gamma_, sigma_, energies_)`</span>"


        :   The energy-resolved, single-shot sibling of [`#!python custom_two()`][calculation-custom_two]: evaluated
            directly by [KITEx][kitex] at fixed broadening and a chosen list of energies, analogous to how
            [`#!python singleshot_conductivity_dc()`][calculation-singleshot_conductivity_dc] relates to
            [`#!python conductivity_dc()`][calculation-conductivity_dc] — no [KITE-tools][kitetools] post-processing step
            is required.

            **Parameters**

            :   | Parameter                                        | Description                                                                                                          |
                |-------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
                | `#!python stream_`:*`#!python list(kite.custom.Vertex, kite.custom.Vertex)`* | The two vertices `#!python [A, B]` defining the operators $A$ and $B$.                            |
                | `#!python num_random_`:*`#!python int`*               | Number of random vectors to use for the stochastic evaluation of the trace.                                          |
                | `#!python num_disorder_`:*`#!python int`*             | Number of different disorder realisations.                                                                           |
                | `#!python gamma_`:*`#!python float`*                  | Imaginary broadening (self-energy) parameter, as in [`#!python singleshot_conductivity_dc()`][calculation-singleshot_conductivity_dc]'s `#!python eta`. |
                | `#!python sigma_`:*`#!python float`*                  | Width of the Gaussian broadening applied to the reconstructed spectral quantity.                                     |
                | `#!python energies_`:*`#!python array_like`*          | List of energies at which `#!python custom_singleshot_two()` will be calculated.                                     |

            !!! Warning "`energies_`/`gamma_`/`sigma_` are not shifted by `energy_shift`"

                Unlike [`#!python singleshot_conductivity_dc()`][calculation-singleshot_conductivity_dc], which subtracts
                [`#!python config.energy_shift`][configuration-energy_shift] from `#!python energy` before exporting it,
                `#!python custom_singleshot_two()`'s `#!python energies_`/`#!python gamma_`/`#!python sigma_` are only
                divided by `#!python config.energy_scale` (this division happens on the KITEx/C++ side, not in Python) —
                `#!python energy_shift` is never subtracted. If your `#!python spectrum_range` is not symmetric about
                zero, `#!python energies_` is measured relative to the Hamiltonian's raw zero, not the center of your
                configured `#!python spectrum_range`.

            !!! Note "No worked example yet"

                There is currently no example script that uses `#!python custom_singleshot_two()`.

## kite.visualize

:   Pure-Python (numpy/matplotlib/scipy) tools for inspecting a native `#!python kite.lattice.Lattice`
    directly — unit cell, Brillouin zone, k-paths, and tight-binding band structure — with **no
    dependency on KITEx/KITE-tools or pybinding**. Meant as a fast sanity check of a lattice definition
    (wrong sign, wrong hopping, wrong lattice vector) before paying for an expensive KPM/Chebyshev run,
    not as a replacement for KITE's own stochastic machinery. Import as `#!python from kite import
    visualize as viz`.

    !!! declaration-function "<span id="visualize-make_path">*function* `#!python viz.make_path(*points, step=0.1, point_labels=None)`</span>"

        :   Build a piecewise-linear k-path through a sequence of high-symmetry points, with point
            density set by each segment's actual Cartesian reciprocal-space distance (not a fixed count
            per segment, and not fractional-coordinate distance). `#!python points` must be given in
            Cartesian reciprocal-space coordinates — the same units `#!python
            lattice.reciprocal_vectors()` returns — not fractional/reduced coordinates (a native
            `#!python kite.lattice.Lattice` method; no separate API doc page exists for this class yet).
            Returns `#!python (k_points, tick_indices, tick_labels)`; `#!python tick_indices`
            locates each input high-symmetry point's row in `#!python k_points`, for tick placement on a
            band-structure plot's x-axis.

    !!! declaration-function "<span id="visualize-plot_unit_cell">*function* `#!python viz.plot_unit_cell(lattice, repeat=(1, 1), ax=None, sublattice_colors=None)`</span>"

        :   Scatter each sublattice's position (stably colored/labeled by name) and draw a line for every
            hopping term to its actual neighbor (including across unit-cell boundaries), repeated over a
            small window of neighboring cells so boundary-crossing bonds are visible. 2D lattices only.

    !!! declaration-function "<span id="visualize-plot_brillouin_zone">*function* `#!python viz.plot_brillouin_zone(lattice, ax=None, k_path=None, length_unit="nm")`</span>"

        :   Draw the 2D Brillouin zone (the Voronoi cell of the reciprocal lattice around $\mathbf k=0$,
            via `#!python lattice.brillouin_zone()`) plus the reciprocal lattice vectors.
            `#!python k_path` optionally overlays a path — either raw points, or
            `#!python make_path(...)`'s full 3-tuple return value directly (in which case high-symmetry
            points get marked and labeled). 2D lattices only.

    !!! declaration-function "<span id="visualize-hamiltonian_k">*function* `#!python viz.hamiltonian_k(lattice, k)`</span>"

        :   Build the Bloch tight-binding Hamiltonian $H(\mathbf k)$ directly from a lattice's own
            hoppings — one sublattice = one orbital (`#!python kite.lattice.Lattice` only stores a scalar
            onsite/hopping value per sublattice pair; multi-orbital-per-sublattice matrices are a separate
            `#!python kite.custom` concern for building observables, not this base Hamiltonian).

            !!! Info "Gauge convention (derived, not assumed)"

                Matches KITE's own ARPES plane-wave state (the full atomic position, cell + sublattice
                offset, enters the phase — see [Spectral Function][spectral-function-example]). For a
                hopping family with `#!python relative_index=R`, `#!python from_sub=a`, `#!python
                to_sub=b`, `#!python energy=t`:
                $$
                H_{ba}(\mathbf k)\mathrel{+}= t\;e^{-i\,\mathbf k\cdot(\mathbf R_{\text{cart}}+\mathbf d_b-\mathbf d_a)}
                $$
                **with a minus sign in the exponent** — get this backwards and centrosymmetric bands
                (e.g. graphene) still look right while anything valley/Rashba/Weyl-asymmetric comes out
                mirror-flipped. Since `#!python add_one_hopping` only ever stores one bond direction, this
                function builds the one-directional sum $H_0(\mathbf k)$ and returns $H_0(\mathbf
                k)+H_0(\mathbf k)^\dagger+\text{onsite}$, reproducing the missing reverse-direction terms
                automatically. The result is checked numerically Hermitian before being returned; a
                non-Hermitian result raises `#!python RuntimeError` rather than being silently returned.

            !!! Example "Verified against several known models, not just graphene"

                `#!python examples/band_structure_graphene.py` checks the three textbook graphene
                energies directly: $E(\Gamma)=\pm3t$, $E(M)=\pm t$, $E(K)=0$ (the Dirac point, exactly
                two-fold degenerate) — all three reproduced to floating-point precision. Also checked
                (not shipped as example scripts, but confirmed while building this feature): the Haldane
                model's gap ($2\times3\sqrt3\,t_2\sin\phi$) matches exactly at both inequivalent BZ
                corners; `#!python kite.repository.group6_tmd.monolayer_3band('MoS2')` (3 coincident-
                position sublattices) gives a clean ~1.85 eV direct gap at K with no accidental
                degeneracy; the T-symmetric cubic Weyl semimetal (`#!python examples/weyl_lt.py`, 3D, two
                sublattices, hoppings given by the model with same-$\mathbf R$ forward/backward terms
                rather than relying on the automatic Hermitian-conjugate addition) reproduces its Weyl
                node exactly at $(\pi/2,\pi/2,\pi/2)$, matching the closed-form
                $H(\mathbf k)=t[\cos k_x\,\sigma_x+\cos k_y\,\sigma_y+\cos k_z\,\sigma_z]$ to machine
                precision — confirming the automatic $H_0+H_0^\dagger$ construction is safe even when a
                lattice's own hopping list already contains what looks like a manually-added return path.
                `#!python kite.repository.phosphorene.monolayer_4band()` (4 sublattices with an
                out-of-plane buckling offset) gives a physically sane, anisotropic, direct ~1.5 eV gap at
                $\Gamma$ — this case caught a real bug (a shape-mismatch crash whenever a sublattice
                position has more components than the lattice has primitive vectors), now fixed by
                zero-padding the in-plane cell-translation vector and truncating the dot product with
                $\mathbf k$ to $\mathbf k$'s own dimensionality, exactly mirroring
                [`#!python plot_unit_cell`][visualize-plot_unit_cell]'s top-down-projection convention for
                buckled lattices.

    !!! declaration-function "<span id="visualize-compute_bands">*function* `#!python viz.compute_bands(lattice, k_path)`</span>"

        :   Diagonalize [`#!python hamiltonian_k`][visualize-hamiltonian_k] at every point of a k-path
            (raw points, or `#!python make_path(...)`'s 3-tuple return value). Returns eigenvalues, shape
            `#!python (n_k, lattice.nsub)`, ascending order.

    !!! declaration-function "<span id="visualize-plot_bands">*function* `#!python viz.plot_bands(lattice, k_path, ax=None, ylabel="Energy (eV)")`</span>"

        :   Plot the band structure along a k-path. If `#!python make_path(...)`'s 3-tuple is passed,
            high-symmetry points get vertical grid lines and x-axis tick labels.

## make_pybinding_model

:   !!! declaration-function "*function* `#!python kite.make_pybinding_model(lattice, disorder=None, disorder_structural=None, shape=None)`"
            
            
    :   Build a Pybinding model with disorder used in Kite. Bond disorder or magnetic field are not currently supported.

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |
            | `#!python shape`:*[`#!python pb.Shape`][shape]*                                            | Optional argument [`#!python pb.Shape`][shape].                                                                                                                                                            |

## estimate_bounds

:   !!! declaration-function "*function* `#!python kite.estimate_bounds(lattice, disorder=None, disorder_structural=None)`"
            
            
    :   Estimate the bounds for the energy, given the [`#!python pb.Lattice`][lattice] and/or [`#!python kite.Disorder`][disorder] and/or [`#!python kite.StructuralDisorder`][structural_disorder].

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |

## config_system

:   !!! declaration-function "*function* `#!python kite.config_system(lattice, config, calculation, modification=None, filename="kite_config.h5", disorder=None, disorder_structural=None)`"
            
            
    :   Export the lattice and related parameters to the *.h5 file.

        **Parameters**
    
        :   | Parameter                                                                                  | Description                                                                                                                                                                                                |
            |--------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
            | `#!python lattice`:*[`#!python pb.Lattice`][lattice]*                                      | Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and onsite energies.                                                                            |
            | `#!python config`:*[`#!python kite.Configuration`][configuration]*                         | [`#!python kite.Configuration`][configuration] object, basic parameters defining size, precision, energy scale and number of decomposition parts in the calculation.                                       |
            | `#!python calculation`:*[`#!python kite.Calculation`][calculation]*                        | [`#!python kite.Calculation`][calculation] object that defines the requested functions for the calculation.                                                                                                |
            | `#!python modification`:*[`#!python kite.Modification`][modification]*                     | If specified [`#!python kite.Modification`][modification] object, has the magnetic field selection, either in terms of field, or in the number of flux quantum through the selected system.                |
            | `#!python filename`: *`#!python str`*                                                      | Filename for the output HDF5-file.                                                                                                                                                                         |
            | `#!python disorder`:*[`#!python kite.Disorder`][disorder]*                                 | Class that introduces [`#!python kite.Disorder`][disorder] into the initially built lattice. For more info check the [`#!python kite.Disorder`][disorder] class.                                           |
            | `#!python disorder_structural`:*[`#!python kite.StructuralDisorder`][structural_disorder]* | Class that introduces [`#!python kite.StructuralDisorder`][structural_disorder] into the initially built lattice. For more info check the [`#!python kite.StructuralDisorder`][structural_disorder] class. |

## LoudDeprecationWarning

:   Deprecationwarning.


[comment]: <> (Classes from Pybinding)
[lattice]: https://docs.pybinding.site/en/stable/_api/pybinding.Lattice.html
[pybinding]: https://docs.pybinding.site/
[shape]: https://docs.pybinding.site/en/stable/api.html#shapes

[comment]: <> (Class StructuralDisorder)
[structural_disorder]: #structuraldisorder
[comment]: <> (Class Parameters)
[structuraldisorder-lattice]: #structuraldisorder-lattice
[structuraldisorder-concentration]: #structuraldisorder-concentration
[structuraldisorder-position]: #structuraldisorder-position
[comment]: <> (Class Methods)
[structuraldisorder-add_vacancy]: #structuraldisorder-add_vacancy
[structuraldisorder-add_structural_disorder]: #structuraldisorder-add_structural_disorder
[structuraldisorder-add_local_vacancy_disorder]: #structuraldisorder-add_local_vacancy_disorder
[structuraldisorder-add_local_bond_disorder]: #structuraldisorder-add_local_bond_disorder
[structuraldisorder-add_local_onsite_disorder]: #structuraldisorder-add_local_onsite_disorder
[structuraldisorder-map_the_orbital]: #structuraldisorder-map_the_orbital

[comment]: <> (Class Disorder)
[disorder]: #disorder
[comment]: <> (Class Parameters)
[disorder-lattice]: #disorder-lattice
[comment]: <> (Class Methods)
[disorder-add_disorder]: #disorder-add_disorder
[disorder-add_local_disorder]: #disorder-add_local_disorder

[comment]: <> (Class Modification)
[modification]: #modification
[comment]: <> (Class Parameters)
[modification-par-magnetic_field]: #modification-par-magnetic_field
[modification-par-flux]: #modification-par-flux
[comment]: <> (Class Attributes)
[modification-atr-magnetic_field]: #modification-atr-magnetic_field
[modification-atr-flux]: #modification-atr-flux

[comment]: <> (Class Calculation)
[calculation]: #calculation
[comment]: <> (Class Parameters)
[calculation-configuration]: #calculation-configuration
[comment]: <> (Class Attributes)
[calculation-get_dos]: #calculation-get_dos
[calculation-get_ldos]: #calculation-get_ldos
[calculation-get_ldos_map]: #calculation-get_ldos_map
[calculation-get_spectral_map]: #calculation-get_spectral_map
[calculation-get_arpes]: #calculation-get_arpes
[calculation-get_gaussian_wave_packet]: #calculation-get_gaussian_wave_packet
[calculation-get_localized_wave_packet]: #calculation-get_localized_wave_packet
[calculation-get_conductivity_dc]: #calculation-get_conductivity_dc
[calculation-get_conductivity_optical]: #calculation-get_conductivity_optical
[calculation-get_conductivity_optical_nonlinear]: #calculation-get_conductivity_optical_nonlinear
[calculation-get_singleshot_conductivity_dc]: #calculation-get_singleshot_conductivity_dc
[calculation-get_custom_one]: #calculation-get_custom_one
[calculation-get_custom_one_local]: #calculation-get_custom_one_local
[calculation-get_custom_two]: #calculation-get_custom_two
[calculation-get_custom_two_local]: #calculation-get_custom_two_local
[calculation-get_custom_ss_two]: #calculation-get_custom_ss_two
[comment]: <> (Class Methods)
[calculation-dos]: #calculation-dos
[calculation-ldos]: #calculation-ldos
[calculation-ldos_map]: #calculation-ldos_map
[calculation-spectral_map]: #calculation-spectral_map
[calculation-arpes]: #calculation-arpes
[calculation-gaussian_wave_packet]: #calculation-gaussian_wave_packet
[calculation-localized_wave_packet]: #calculation-localized_wave_packet
[calculation-conductivity_dc]: #calculation-conductivity_dc
[calculation-conductivity_optical]: #calculation-conductivity_optical
[calculation-conductivity_optical_nonlinear]: #calculation-conductivity_optical_nonlinear
[calculation-singleshot_conductivity_dc]: #calculation-singleshot_conductivity_dc
[calculation-add_orbital_index]: #calculation-add_orbital_index
[calculation-add_orbital_coupling]: #calculation-add_orbital_coupling
[calculation-custom_one]: #calculation-custom_one
[calculation-custom_one_local]: #calculation-custom_one_local
[calculation-custom_two]: #calculation-custom_two
[calculation-custom_two_local]: #calculation-custom_two_local
[calculation-custom_singleshot_two]: #calculation-custom_singleshot_two

[comment]: <> (Class Configuration)
[configuration]: #configuration
[comment]: <> (Class Parameters)
[configuration-divisions]: #configuration-divisions
[configuration-length]: #configuration-length
[configuration-boundaries]: #configuration-boundaries
[settings-boundaries]: ../documentation/settings.md#boundaries
[configuration-is_complex]: #configuration-is_complex
[configuration-precision]: #configuration-precision
[configuration-spectrum_range]: #configuration-spectrum_range
[configuration-angles]: #configuration-angles
[configuration-custom_local]: #configuration-custom_local
[configuration-custom_local_print]: #configuration-custom_local_print
[comment]: <> (Class Attributes)
[configuration-energy_scale]: #configuration-energy_scale
[configuration-energy_shift]: #configuration-energy_shift
[configuration-comp]: #configuration-comp
[configuration-prec]: #configuration-prec
[configuration-div]: #configuration-div
[configuration-bound]: #configuration-bound
[configuration-leng]: #configuration-leng
[configuration-type]: #configuration-type
[configuration-custom_pot]: #configuration-custom_pot
[configuration-print_custom_pot]: #configuration-print_custom_pot
[comment]: <> (Class Methods)
[configuration-set_type]: #configuration-set_type

[make_pybinding_model]: #make_pybinding_model
[estimate_bounds]: #estimate_bounds
[config_system]: #config_system
[loud_deprecation_warning]: #louddeprecationwarning

[kitex]: kitex.md
[kitetools]: kite-tools.md
[kite-tools-customone]: kite-tools.md#kite-tools-customone
[kitetools-output]: kite-tools.md#output
[HDF5]: https://www.hdfgroup.org
[tutorial-hdf5]: ../documentation/editing_hdf_files.md

[magnetic-field]: ../documentation/magnetic.md
[custom-vertex-example]: ../documentation/examples/custom_vertex_operators.md
[rashba-edelstein-example]: ../documentation/examples/rashba_edelstein.md
[spectral-function-example]: ../documentation/examples/spectral_function.md

[^1]: A. Weiße, G. Wellein, A. Alvermann, and H. Fehske, [Rev. Mod. Phys. 78, 275 (2006)](https://doi.org/10.1103/RevModPhys.78.275).
[^2]: S. M. João, M. Anđelković, L. Covaci, T. G. Rappoport, João M. Viana Parente Lopes, and A. Ferreira, [R. Soc. open sci. 7, 191809 (2020)](https://royalsocietypublishing.org/doi/10.1098/rsos.191809).
[^3]: H. P. Veiga, D. R. Pinheiro, J. P. Santos Pires, and J. M. Viana Parente Lopes, "Markov Inequality as a Tool for Linear-Scaling Estimation of Local Observables," [Phys. Rev. Research, doi 10.1103/qb1w-44r1](https://journals.aps.org/prresearch/abstract/10.1103/qb1w-44r1), [arXiv:2510.21688](https://arxiv.org/abs/2510.21688).