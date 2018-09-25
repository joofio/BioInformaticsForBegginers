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

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)


def Pr(Text, Profile):
    k = len(Profile['A'])
    result = 1
    i = 0
    for symbol in Text:
        result = result*Profile[symbol][i]
        i += 1
    return result

# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats


def ProfileMostProbableKmer(text, k, profile):
    n = len(text)
    result = text[0:k]
    m = 0
    for i in range(n-k+1):
        if Pr(text[i:i+k], profile) > m:
            result = text[i:i+k]
            m = Pr(text[i:i+k], profile)
    return result

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
