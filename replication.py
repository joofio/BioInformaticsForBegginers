##Find a Pattern and the number of times it appears
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

# Copy your FrequencyMap() function here.


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq.keys():
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq

# Input:  A string Pattern
# Output: The reverse of Pattern


def Reverse(Pattern):
    rev = ''
    for char in reversed(Pattern):
        rev += char
    return rev
    # your code here

# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).


def Complement(Pattern):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    comp = ''
    for char in Pattern:
        comp += complement[char]
    return comp


def reversecompl(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

#Return the positions of the Pattern inside a string of genome


def PatternMatching(Pattern, Genome):
    positions = []  # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions


#Count of symbol in the half genome starting at position X
def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
        array.append(count)
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    skew=[]
    skew.append(0)
    for i in range(len(Genome)):
        if  Genome[i] =='G':
            skew.append(skew[i]+1)
        elif  Genome[i] =='C':
            skew.append(skew[i]-1)
        else:
            skew.append(skew[i])
    return skew
         




