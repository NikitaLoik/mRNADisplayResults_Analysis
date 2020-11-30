#String = 'ab(cd)2e((fg)2h(ij))(klmnopq)2rstuv'
String = 'q(abc)'

# define ExpandedString
ExpandedString = String
# define InnerMostSubString = []
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

print ExpandedString
