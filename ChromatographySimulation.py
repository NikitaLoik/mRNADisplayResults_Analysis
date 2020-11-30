import matplotlib.pyplot as plt

#   store information about Partitioning Coefficient (PC) and Analyte Quantity
class Analyte:
    AnalyteCounter = 0
    
    def __init__(self, Name, PC, AnalyteQuantity):
        self.Name = Name
        self.PC = PC
        self.AnalyteQuantity = AnalyteQuantity
        Analyte.AnalyteCounter += 1
        
    def StPhaseFraction(self):
        StPhaseFraction = 1/(1 + self.PC)
        return StPhaseFraction
    def MPhaseFraction(self):
        MPhaseFraction = self.PC/(1 + self.PC)
        return MPhaseFraction
    def ShowAnalyteCounter(self):
        print "Number of analytes %d" % Analyte.AnalyteCounter
                
        
#    tests    
Proline = Analyte('Proline', 0.5, 1.0)
#    print Proline.Name
#    print Proline.PC
#    print Proline.AnalyteQuantity
#    print Proline.StPhaseFraction()
#    print Proline.MPhaseFraction()
#    print Proline.StPhaseFraction() + Proline.MPhaseFraction()

#    print Proline.ShowAnalyteCounter()

Lysine = Analyte('Lysine', 1, 1.0)
#    print Lysine.Name
#    print Lysine.PC
#    print Lysine.AnalyteQuantity
#    print Lysine.StPhaseFraction()
#    print Lysine.MPhaseFraction()
#    print Lysine.StPhaseFraction() + Lysine.MPhaseFraction()

#    print Proline.ShowAnalyteCounter()
#    print Lysine.ShowAnalyteCounter()


#------------------------------------------------------------------------------#
        
#   store information about the Number of Theoretical Plates (NTP)
class Column:
    def __init__(self, NTP):
        self.NTP = NTP
#    tests    
NewColumn = Column(100.0)
#    print NewColumn.NTP
#------------------------------------------------------------------------------#
#   store information about the Run Time
class Separation:
    def __init__(self, RunTime):
        self.RunTime = RunTime
#    tests    
NewSeparation = Separation(100.0)
#    print NewSeparation.RunTime
#------------------------------------------------------------------------------#

class Chromatogram(Analyte, Column, Separation):
    #   inheritance
    def __init__(self, Name, PC, AnalyteQuantity, NTP, RunTime):
        Analyte.__init__(self, Name, PC, AnalyteQuantity)
        Column.__init__(self, NTP)
        Separation.__init__(self, RunTime)    
    
    def Chromatography(self):

        #   defining zero point on chromatogram
        Chromatogram = [0.0]
        
        #   defining precursory starting time point
        PrecursoryTimePoint = [self.AnalyteQuantity]+[0.0]*(int(self.NTP)-1)
        #   defining current time point
        CurrentTimePoint = [0.0]*int(self.NTP)
                       
        for TimePoint in range(int(self.RunTime)):            
            CurrentTimePoint[0] = PrecursoryTimePoint[0]*self.StPhaseFraction()
            for TP in range(1, int(self.NTP)):
                CurrentTimePoint[TP] = PrecursoryTimePoint[TP-1]*self.MPhaseFraction() + PrecursoryTimePoint[TP]*self.StPhaseFraction()
            PrecursoryTimePoint = CurrentTimePoint[:]                
            Chromatogram = Chromatogram + [CurrentTimePoint[-1]]
        return Chromatogram
        
    def NoisyChromatography(self):
        import random

        #   defining zero point on chromatogram
        Chromatogram = [0.0]
        SNRatio = 3
        NoisyChromatogram = []
        
        #   defining precursory starting time point
        PrecursoryTimePoint = [self.AnalyteQuantity]+[0.0]*(int(self.NTP)-1)
        #   defining current time point
        CurrentTimePoint = [0.0]*int(self.NTP)
                       
        for TimePoint in range(int(self.RunTime)):            
            CurrentTimePoint[0] = PrecursoryTimePoint[0]*self.StPhaseFraction()
            for TP in range(1, int(self.NTP)):
                CurrentTimePoint[TP] = PrecursoryTimePoint[TP-1]*self.MPhaseFraction() + PrecursoryTimePoint[TP]*self.StPhaseFraction()
            PrecursoryTimePoint = CurrentTimePoint[:]                
            Chromatogram = Chromatogram + [CurrentTimePoint[-1]]
            
        for TimePoint in range(len(Chromatogram)):
            NoisyTimePoint = Chromatogram[TimePoint] + max(Chromatogram)*random.random()*(1.0/SNRatio)
            NoisyChromatogram = NoisyChromatogram + [NoisyTimePoint]
                    
        return NoisyChromatogram        
        
    def RetentionTime(self):
        RetentionTime = 0
        Chromatogram = self.Chromatography()
        for TimePoint in range(len(Chromatogram)):
            if Chromatogram[TimePoint] == max(Chromatogram):
                RetentionTime = TimePoint
        return RetentionTime
        
    def NoisyRetentionTime(self):
        NoisyRetentionTime = 0
        NoisyChromatogram = self.NoisyChromatography()
        for TimePoint in range(len(NoisyChromatogram)):
            if NoisyChromatogram[TimePoint] == max(NoisyChromatogram):
                NoisyRetentionTime = TimePoint
        return NoisyRetentionTime       
        
    def PeakWidthAtHalfHeight(self):
        Chromatogram = self.Chromatography()
#        RetentionTime = self.RetentionTime()
#        StartTime = 0
#        EndTime = 0
        PeakHalfHeight = 0.5*max(Chromatogram)
        count = 0
        for TimePoint in range(len(Chromatogram)):
            if (Chromatogram[TimePoint] - PeakHalfHeight) > 0.0:
                count =+ 1
        return count
                   
    def ChromatogramPlot(self):
        plt.plot(self.Chromatography())        
        return plt.show()
        
    def NoisyChromatogramPlot(self):
        plt.plot(self.NoisyChromatography())        
        return plt.show()

#class ExtractedChromatogram():
    


#    tests    
ProlineChromatogram = Chromatogram('Proline', 0.2, 1.0, 100.0, 1000.0)
LysineChromatogram = Chromatogram('Lysine', 0.5, 1.0, 100.0, 1000.0)
print ProlineChromatogram.ChromatogramPlot()
print LysineChromatogram.ChromatogramPlot()
print ProlineChromatogram.RetentionTime()
print LysineChromatogram.RetentionTime()
print ProlineChromatogram.PeakWidthAtHalfHeight()
print LysineChromatogram.PeakWidthAtHalfHeight()

print ProlineChromatogram.NoisyChromatogramPlot()
print LysineChromatogram.NoisyChromatogramPlot()
print ProlineChromatogram.NoisyRetentionTime()
print LysineChromatogram.NoisyRetentionTime()

#------------------------------------------------------------------------------#
