import numpy as np

def global_sequence_alignment(seq_a, seq_b, match=1, mismatch=-1, gap_penalty=-2):
    len_a, len_b = len(seq_a), len(seq_b)
    score_matrix = np.zeros((len_a + 1, len_b + 1))

    for i in range(len_a + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(len_b + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, len_a + 1):
        for j in range(1, len_b + 1):
            if seq_a[i - 1] == seq_b[j - 1]:
                match_score = score_matrix[i - 1][j - 1] + match
            else:
                match_score = score_matrix[i - 1][j - 1] + mismatch

            delete_score = score_matrix[i - 1][j] + gap_penalty
            insert_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(match_score, delete_score, insert_score)

    align_a, align_b = "", ""
    i, j = len_a, len_b
    while i > 0 or j > 0:
        if i > 0 and j > 0 and seq_a[i - 1] == seq_b[j - 1]:
            align_a = seq_a[i - 1] + align_a
            align_b = seq_b[j - 1] + align_b
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            align_a = seq_a[i - 1] + align_a
            align_b = "-" + align_b
            i -= 1
        else:
            align_a = "-" + align_a
            align_b = seq_b[j - 1] + align_b
            j -= 1

    return score_matrix[len_a][len_b], align_a, align_b

sequence_a = input("Enter the first string: ")
sequence_b = input("Enter the second string: ")
score, aligned_a, aligned_b = global_sequence_alignment(sequence_a, sequence_b)
print("Similarity Score:", score)
print("Aligned Sequence A:", aligned_a)
print("Aligned Sequence B:", aligned_b)

