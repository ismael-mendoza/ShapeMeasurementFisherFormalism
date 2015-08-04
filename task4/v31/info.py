import defaults

import os

class Info: 
    """Contains what is written in the information text file and printed out 
    with verbose option.
    """
    
    def __init__(self, g_parameters, fish = None, number_fits = None):
        params = g_parameters.params 
        self.galaxy = []
        self.fisher = []
        self.fits = []

        #contains all the lines of galaxy info to be writen into the file. 
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
            +(
            'the results were produced with:',
            os.path.dirname(os.path.realpath(__file__))
            )
        )

        if fish:
            self.fisher = list(
                (
                '',
                'Fisher analysis (and fittings if included) used the following snr: ' + str(fish.snr),
                '',
                'Steps used for the derivatives: '
                ) 
                + 
                tuple(param + ': ' + str(fish.steps[param]) for param in fish.steps.keys())
            )

        if number_fits:
            self.fits = list(
                (
                'Number of fits: ' + str(number_fits)
                )
            )

        self.lines = self.galaxy + self.fisher + self.fits

    def dictToStringList(self,dct): 
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
                txtfile.write(line + '\n') #last character to skip lines.
