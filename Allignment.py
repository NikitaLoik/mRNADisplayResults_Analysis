StringA = 'abcdefghij'
StringB = 'oiupmcdmw'

print len(StringA)
print len(StringB)

print StringA.find('cd')
print StringB.find('cd')
Correction = StringA.find('cd') - StringB.find('cd') + 10

BlankString = '__________'
TopString = BlankString + StringA

BalancingString = Correction * '_'
String = BalancingString + StringB

print TopString
print String