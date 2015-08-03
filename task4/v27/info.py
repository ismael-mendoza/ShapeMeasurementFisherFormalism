import defaults

import os

class Info: 
    """Contains what is written in the information text file and printed out 
    with verbose option.
    """
    
    def __init__(self, g_parameters, fish = None, fits_biases = None, number_fits = None):
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
            'pixel_scale: ' + str(pixel_scale),
            '',
            'Galaxies drawn have the following parameters:'
            ) 
            + 
            tuple(param + ': ' + str(params[param]) for param in params.keys()
                 ) 
        )

        if fish:
            self.fisher = list(
                (
                '',
                'Fisher matrix elements:',
                )
                + 
                tuple(dictToStringList(fish.fisher_matrix))
                + (
                '',
                'Covariance matrix elements:'
                )
                + 
                tuple(dictToStringList(fish.covariance_matrix))
                +
                (
                '',
                'Biases for each parameter'
                )
                + 
                tuple(dictToStringList(fish.biases))
                + 
                (
                '',
                'Fisher analysis used the following snr): ' + str(fish.snr),
                '',
                'Steps used for the derivatives: '
                ) 
                + 
                tuple(param + ': ' + str(fish.steps[param]) for param in fish.steps.keys())
            )

        #add covariances in the future? (unlikely)
        if fits_biases and number_fits:
            self.fits = list(
                (
                'Number of fits: ' + str(number_fits),
                'Biases were obtained by averaging over fits of noisy instantiations of the same image: '
                )
                +
                tuple(dictToStringList(fits_biases))
            )

        self.lines = self.galaxy + self.fisher + self.fits

    def dictToStringList(dct): 
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

    def writeInfo(self, project_dir, info_file):
 
        with open(os.path.join(project_dir, info_file), 'w') as txtfile:
            for line in self.galaxy + self.fisher:
                txtfile.write(line + '\n') #last character to skip lines.
