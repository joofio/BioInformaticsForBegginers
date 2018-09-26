# Input:  A set of kmers Motifs
# Output: Count(Motifs)
import random

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
    result = 1
    i = 0
    for symbol in Text:
        result = result*Profile[symbol][i]
        i += 1
    return result


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


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}  # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)+4
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j] = profile[symbol][j]/(t)
    return profile    

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    n=len(Profile['A'])
    result=[]
    for item in Dna:
        result.append(ProfileMostProbableKmer(item,n,Profile))
    return result

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    result=[]
    for i in range(t):   
        rnd=random.randint(0, len(Dna[i])-k)
        result.append(Dna[i][rnd:rnd+k])
    return result

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    sum=0
    for i in Probabilities.values():
        sum=sum+i
    for symbol in Probabilities.keys():
        Probabilities[symbol]=Probabilities[symbol]/sum
    return Probabilities

# first, import the random package
import random
# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    prob=0
    rnd=random.uniform(0, 1)
    for i in Probabilities:
       # print(str(prob) +' '+str(prob+Probabilities[i]))
        if prob<rnd<prob+Probabilities[i]:
            kmer=i
            break
        else:
            prob=prob+Probabilities[i]
    return kmer

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
# TODO
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k]) 
    
    for j in range(1,N):
        i=random.randint(0,t)
        n = len(Dna[0])
        for i in range(n-k+1):
            Motifs = []
            Motifs.append(Dna[0][i:i+k])
            for j in range(1, t):
                P = ProfileWithPseudocounts(Motifs[0:j])
                Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
        return BestMotifs