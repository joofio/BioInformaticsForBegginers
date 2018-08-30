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
