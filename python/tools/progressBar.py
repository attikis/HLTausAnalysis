###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################
###############################################################
### All imported modules
###############################################################
### System modules
import progressbar

###############################################################
### Define class here
###############################################################
class ProgressBar(object): 
    def __init__(self, maxValue):
        self.widgets  = [progressbar.FormatLabel(''), ' ', progressbar.Percentage(), ' ', progressbar.Bar('/'), ' ', progressbar.RotatingMarker()]
        self.pBar     = progressbar.ProgressBar(widgets=self.widgets, maxval=maxValue)
        self.Start()

    def Start(self):
        '''
        Simple module to create and initialise a progress bar. The argument "maxvalue" refers to the 
        total number of tasks to be completed. This must be defined at the start of the progress bar.
        '''
        #if self.pBar.start_time is None:
        self.pBar.start()
        return 

    def Update(self, counter):
        '''
        Update the progress bar.
        '''
        self.pBar.update(counter)
        return 

    def Finish(self):
        '''
        Finish the progress bar.
        '''
        self.pBar.finish()
        return 
