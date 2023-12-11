import numpy as np

def viterbi_alignment(seq1, seq2, transition_prob, emission_prob):
    len1, len2 = len(seq1), len(seq2)
    states = ['MATCH', 'INSERT', 'DELETE']
    num_states = len(states)

    dp = np.full((len1 + 1, len2 + 1, num_states), float('-inf'))
    path = np.zeros((len1 + 1, len2 + 1, num_states), dtype=int)

    dp[0][0][0] = 1.0

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            for k in range(num_states):
                for l in range(num_states):
                    prob = dp[i-1][j-1][l] * transition_prob[l][k]

                    if k == 0:
                        if seq1[i-1] == seq2[j-1]:
                            prob *= emission_prob[k][symbol_to_index(seq1[i-1])]
                        else:
                            prob = float('-inf')
                    else:
                        prob *= emission_prob[k][symbol_to_index(seq1[i-1])]

                    if prob > dp[i][j][k]:
                        dp[i][j][k] = prob
                        path[i][j][k] = l

    aligned_seq1 = ""
    aligned_seq2 = ""
    state_sequence = []
    i, j, state = len1, len2, np.argmax(dp[len1][len2])

    while i > 0 and j > 0:
        state_sequence.append(states[state])
        if state == 0:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif state == 1:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        state = path[i][j][state]

    return list(reversed(state_sequence)), seq1, seq2, aligned_seq1, aligned_seq2


def symbol_to_index(symbol):
    return {'A': 0, 'D': 1}.get(symbol, -1)

transition_prob = [
    [0.9, 0.05, 0.05],
    [0.45, 0.1, 0.45],
    [0.45, 0.45, 0.1]
]

emission_prob = [
    [0.6, 0.4],
    [0.5, 0.5],
    [0.5, 0.5]
]

seq1 = "ADADA"
seq2 = "AADDD"

viterbi_states, original_seq1, original_seq2, aligned_seq1, aligned_seq2 = viterbi_alignment(seq1, seq2, transition_prob, emission_prob)

print("States: ")
print(viterbi_states)
print("Original Sequences: ")
print("S1: ", original_seq1)
print("S2: ", original_seq2)
print("Alignen Sequences: ")
print("S1: ", aligned_seq1)
print("S2: ", aligned_seq2)
