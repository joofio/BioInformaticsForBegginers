# Input:  A set of kmers Motifs
# Output: Count(Motifs)


def Count(Motifs):
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j] = profile[symbol][j]/(t)
    return profile


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.


def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    count = 0
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] != consensus[j]:
                count += 1
    return count