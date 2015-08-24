import defaults

import os

class Info:
    """Contains what is written in the information text file and printed out
    with verbose option in the interface files.


    The :class:`Info` object will contain different information depending on
    the different parameters that are passed to it. All of its attributes are
    list of strings that are printed in a format specified when the object is
    created.

        Args:
            g_parameters(:class:`GParameters`): Object that contains the
                                                relevant information of the
                                                parameters of the galaxies
                                                used.
            fish(:class:`Fisher`): Object containing the results of the
                                   different fisher analysis done over the
                                   specified galaxies.
            number_fits(int): Number of fits done if a


        Attributes:


        Returns:
            A galsim.Image object.
    """

    def __init__(self, g_parameters, fish=None, number_fits=None,
                 init_values=None, minimums=None, maximums=None):
        params = g_parameters.params
        self.galaxy = []
        self.fisher = []
        self.fits = []
        self.init_values = []
        self.minimums = []
        self.maximums = []

        # contains all the lines of galaxy info to be writen into the file.
        self.galaxy = list(
            (
                'Default values used in the analysis:',
                'nx: ' + str(defaults.NX),
                'ny: ' + str(defaults.NY),
                'pixel_scale: ' + str(defaults.PIXEL_SCALE),
                '',
                'Galaxies drawn have the following parameters:'
            )
            +
            tuple(param + ': ' + str(params[param]) for param in params.keys()
                  )
            + (
                'the results were produced with:',
                os.path.dirname(os.path.realpath(__file__))
            )
        )

        if fish:
            self.fisher = list(
                (
                    '',
                    'Fisher analysis (and fittings if included) used the following snr: ' + str(
                        fish.snr),
                    '',
                    'Steps used for the derivatives: '
                )
                +
                tuple(param + ': ' + str(fish.steps[param])
                      for param in fish.steps)
            )

        if number_fits:
            self.fits = list(
                (
                    'Number of fits:',
                    str(number_fits)
                )
            )

        if init_values:
            self.init_values = list(
                ('',
                 'We used the following initial values for the fitting:',)
                +
                tuple(param + ': ' + str(init_values[param]) for param in init_values
                      )

            )

        if minimums:
            self.minimums = list(
                ('',
                 'We used the following minimum values for the fitting:',)
                +
                tuple(param + ': ' + str(minimums[param]) for param in minimums
                      )
            )

        if maximums:
            self.maximums = list(
                ('',
                 'We used the following maximum values for the fitting:',)
                +
                tuple(param + ': ' + str(maximums[param]) for param in maximums
                      )
            )

        self.lines = (self.galaxy + self.fisher + self.fits + self.init_values
                      + self.minimums + self.maximums)

    def dictToStringList(self, dct):
        """Returns a list of strings consisting of the form 'key:value' for a
        given dictionary.
        """
        if type(dct.keys()[0]) == str:
            return [str(key) + ': ' + str(dct[key]) for key in dct.keys()]
        else:
            return [str(','.join(key)) + ': ' + str(dct[key])
                    for key in dct.keys()
                    ]

    def printInfo(self):
        for line in self.lines:
            print line

    def writeInfo(self, project):
        info_file_name = os.path.join(project, defaults.INFO_FILE)
        with open(info_file_name, 'w') as txtfile:
            for line in self.lines:
                txtfile.write(line + '\n')  # last character to skip lines.
