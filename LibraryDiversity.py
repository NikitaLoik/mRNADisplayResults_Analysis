ComplementaryDNA = {'A':'T','C':'G','G':'C','T':'A'}
ExtendedDNA = {'A':'A','C':'C','G':'G','T':'T','R':('A','G'),'Y':('C','T'),'S':('G','C'),'W':('A','T'),'K':('G','T'),'M':('A','C'),'B':('C','G','T'),'D':('A','G','T'),'H':('A','C','T'),'V':('A','C','G'),'N':('A','C','G','T')}
transcription = {'A':'A','C':'C','G':'G','T':'U','U':'T'}
RNA = {'A':'U','C':'G','G':'C','U':'A'}
translation = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N','ACA':'T','ACC':'T','ACG':'T','ACU':'T','AGA':'R','AGC':'S','AGG':'R','AGU':'S','AUA':'I','AUC':'I','AUG':'M','AUU':'I','CAA':'Q','CAC':'H','CAG':'Q','CAU':'H','CCA':'P','CCC':'P','CCG':'P','CCU':'P','CGA':'R','CGC':'R','CGG':'R','CGU':'R','CUA':'L','CUC':'L','CUG':'L','CUU':'L','GAA':'E','GAC':'D','GAG':'E','GAU':'D','GCA':'A','GCC':'A','GCG':'A','GCU':'A','GGA':'G','GGC':'G','GGG':'G','GGU':'G','GUA':'V','GUC':'V','GUG':'V','GUU':'V','UAA':'ochre','UAC':'Y','UAG':'amber','UAU':'Y','UCA':'S','UCC':'S','UCG':'S','UCU':'S','UGA':'opal','UGC':'C','UGG':'W','UGU':'C','UUA':'L','UUC':'F','UUG':'L','UUU':'F'
}

def BracketOpener(String):

    ExpandedString = String
    InnerMostSubString = ''
    ExpandedBracket = ''
    Number = 1
    DigitNumber = 0
    DigitOne = ''
    DigitTwo = ''
    
    while ')' in ExpandedString:
        # define ClosingBracketPosition
        ClosingBracketPosition = ExpandedString.find(')')
        
        if len(ExpandedString) >= ClosingBracketPosition + 3:
            DigitOne = ExpandedString[ClosingBracketPosition + 1]
            DigitTwo = ExpandedString[ClosingBracketPosition + 2]
        elif len(ExpandedString) == ClosingBracketPosition + 2:
            DigitOne = ExpandedString[ClosingBracketPosition + 1]
            DigitTwo = ''
        elif len(ExpandedString) == ClosingBracketPosition + 1:
            DigitOne = ''
            DigitTwo = ''
              
        # define SubString, equal to String beginning after the first opening bracket
        if DigitOne.isdigit() and DigitTwo.isdigit():
            Number = int(DigitOne + DigitTwo)
            DigitNumber = 2
        elif DigitOne.isdigit() and not DigitTwo.isdigit():
            Number = int(DigitOne)
            DigitNumber = 1
        elif not DigitOne.isdigit() and not DigitTwo.isdigit():
            Number = 1
            DigitNumber = 0
        
        SubString = ExpandedString[:ClosingBracketPosition]
        # define OpeningBracketPosition
        OpeningBracketPosition = SubString.rfind('(')
        # define InnerMostSubString = ''
        InnerMostSubString = SubString[OpeningBracketPosition + 1:]
        ExpandedBracket = InnerMostSubString * Number
        ExpandedString = ExpandedString[:OpeningBracketPosition] + ExpandedBracket + ExpandedString[ClosingBracketPosition + 1 + DigitNumber:]
    
    return ExpandedString

SubLibrariesNumber = int(raw_input('Number of SubLibraries '))
if SubLibrariesNumber is not int:
    print 'Are you fucking with me?'

LibraryDiversity = 0

for i in range(SubLibrariesNumber):
    SubLibrary = raw_input('paste the DNA or RNA sequence of the sublibrary ' + str(i+1) + ' ')
    SubLibrary = BracketOpener(SubLibrary)
    SubLibraryCAP = SubLibrary.upper()
    print SubLibraryCAP
    SubLibraryDiversity = 1
    for j in range(len(SubLibraryCAP)):
        if SubLibraryCAP[j] in ExtendedDNA:
            DNAChar = SubLibraryCAP[j]
            SubLibraryDiversity = SubLibraryDiversity * len(ExtendedDNA[DNAChar])
            print DNAChar
            print SubLibraryDiversity
        else:
            DNAChar = transcription[SubLibraryCAP[j]]
            SubLibraryDiversity = SubLibraryDiversity * len(ExtendedDNA[DNAChar])
            print DNAChar
            print SubLibraryDiversity
        
    LibraryDiversity = LibraryDiversity + SubLibraryDiversity
    print 'SubLibrary ' + str(i+1) + ' diversity is ' "%.2g" %(SubLibraryDiversity)
    
print 'Library diversity is ' + "%.2g" %(LibraryDiversity)