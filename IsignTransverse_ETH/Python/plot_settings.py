
from distutils.command.sdist import sdist
import numpy as np

options = ['L', 'J', 'g', 'h', 'w', 'k']
#--- class for settings to control
class plot_settings_class:
    """
    General class to describe behavior of plot structure with additional rescaling options

    Attributes:
    --------------------
    --------------------
    settings:      main class attribute, stores dictionairy with user settings such as:

            * 'vs'              -   value on x-axis, one of the following ['L', 'J', 'g', 'h', 'w', 'k']

            * 'scaling'         -   value changed for each curve (on legend), same possibilities as above

            * 'x_scale'         -   scaling on x-axis: linear, log

            * 'y_scale'         -   scaling on y-axis: linear, log

            * 'physical_units'  -   rescale qunatity to physical units

            * 'rescaleY'        -   boolean value whether to rescale y-values by function given by 'func_y'

            * 'func_y'          -   function type of rescaling; options: 'exp', 'log', 'power-law'; add '_inv' to name to have the inverse, any other option results in no scaling

            * 'nu_y'            -   exponent in function func_y: 'exp'=exp(y^nu_y),     'log'=log_{nu_y}(y),    'power-law'=y^{nu_y}

            * 'rescaleX'        -   same as 'rescaleY' but for x-axis

            * 'func_x'          -   same as 'func_y' but for x-axis

            * 'nu_x'            -   same as 'nu_y'

        #-- automatically set variables
            * 'vs_idx', 'scaling_idx'       -   internal use

            * 'func_x_name', 'func_y_name'  -   label for functions on rescale xy
    """
    settings = {}

    def __init__(self, settings_dict_in):
        self.settings = settings_dict_in
        self.settings['vs_idx'] = options.index(self.settings['vs'])
        self.settings['scaling_idx'] = options.index(self.settings['scaling'])
        #if self.settings['rescaleX'] : self.settings['x_scale'] = 'log'
        
        self.__set_func_names()
    

    #-- PRIVATE METHOD TO SET FUNC NAMES
    def __set_func_names(self):
        dumm = ['x', 'y']
        for ss in dumm:
            is_int = (self.settings['nu_'+ss] - int(self.settings['nu_'+ss])) == 0
            if self.settings['rescale' + ('X' if ss == 'x' else 'Y')]:
                if self.settings['func_' + ss] == 'exp':           self.settings['func_' + ss + '_name'] = r"$exp(Q^{%d})$"%self.settings['nu_'+ss] if is_int else r"$exp(Q^{%.2f})$"%self.settings['nu_'+ss]
                elif self.settings['func_' + ss] == 'exp_inv':     self.settings['func_' + ss + '_name'] = r"$exp(-Q^{%d})$"%self.settings['nu_'+ss] if is_int else r"$exp(-Q^{%.2f})$"%self.settings['nu_'+ss]
                elif self.settings['func_' + ss] == 'log':         self.settings['func_' + ss + '_name'] = r"$log_{%d}\ Q$"%self.settings['nu_'+ss] if is_int else r"$log_{%.2f}\ Q$"%self.settings['nu_'+ss]
                elif self.settings['func_' + ss] == 'log_inv':     self.settings['func_' + ss + '_name'] = r"$log^{-1}_{%d}\ Q$"%self.settings['nu_'+ss] if is_int else r"$log^{-1}_{%.2f}\ Q$"%self.settings['nu_'+ss]
                elif self.settings['func_' + ss] == 'power-law':   self.settings['func_' + ss + '_name'] = r"$Q^{%d}$"%self.settings['nu_'+ss] if is_int else r"$Q^{%.2f}$"%self.settings['nu_'+ss]
                else: self.settings['func_' + ss + '_name'] = r"$Q$"
            else : 
                self.settings['func_' + ss + '_name'] = r"$Q$"

#------ ATTRIBUTES
    @property
    def settings(self):
        return self._settings
    
    @settings.setter
    def settings(self, input_settings):
        self._settings = input_settings

#------ RESCALING FUNCTIONS
    def rescale(self, value, ss = 'x'):
        if self.settings['rescale' + ('X' if ss == 'x' else 'Y')]:
            if self._settings['func_' + ss] == 'exp':           return np.exp( value**self._settings['nu_' + ss])
            elif self._settings['func_' + ss] == 'exp_inv':     return np.exp(-value**self._settings['nu_' + ss])
            elif self._settings['func_' + ss] == 'log':         return np.log(value) / np.log(self._settings['nu_' + ss])
            elif self._settings['func_' + ss] == 'log_inv':     return np.log(self._settings['nu_' + ss]) / np.log(value)
            elif self._settings['func_' + ss] == 'power-law':   return value**self._settings['nu_' + ss]
            else: return value
        else: 
            return value


#------ NON-GENERIC SETTERS
    def set_vs(self, vs_option):
        """
        Set parameter on x-axis, i.e. plot y-data vs ...
        """
        self.settings['vs'] = vs_option
        self.settings['vs_idx'] = options.index(self.settings['vs'])


    def set_scaling(self, scaling_option):
        """
        Set parameter changed for every new plot, i.e. parameter in legend
        """
        self.settings['scaling'] = scaling_option
        self.settings['scaling_idx'] = options.index(self.settings['scaling'])


    def set_scales(self, xscale = 'log', yscale = 'log'):
        """
        Set scales of xy axis: linear or log
        """
        self.settings['x_scale'] = xscale
        self.settings['y_scale'] = yscale


    def set_x_rescale(self, rescale = 0, function = 'power-law', exponent = 2.0):
        """
        Set rescale functions on x-axis
        """
        self.settings['rescaleX'] = rescale
        self.settings['nu_x'] = exponent
        self.settings['func_x'] = function
        self.__set_func_names()


    def set_y_rescale(self, rescale = 0, function = 'power-law', exponent = 2.0):
        """
        Set rescale functions on y-axis
        """
        self.settings['rescaleY'] = rescale
        self.settings['nu_y'] = exponent
        self.settings['func_y'] = function
        self.__set_func_names()