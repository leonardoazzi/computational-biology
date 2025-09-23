# Desenvolva um programa que implemente o algoritmo de Needleman-Wunsch. 
# A implementação deve incluir a criação da matriz de pontuação e a função
# de traceback para reconstruir o alinhamento.

import numpy as np
from enum import Enum
import sys

# Valor para GAP: -3
# Valor para Match: +2
# Valor para Mismatch: -2
class Score(Enum):
    GAP = -3
    MATCH = +2
    MISMATCH = -2
    # MATCH = +5
    # MISMATCH = -3
    # GAP = -4

# Enums de movimentos    
UP = 0
LEFT = 1
DIAG = 2

def read_sequence(tax_id: str) -> str:
    with open(f"./sequences/{tax_id}.txt", "r") as file:
        return file.read().strip()

def get_score(char_1: str, char_2: str) -> int:
    score = 0
    if char_1 == char_2:
        score = Score.MATCH.value
    else:
        score = Score.MISMATCH.value
    
    return score

def smith_waterman(
    seq_1: str,
    seq_2: str,
    start_pos: tuple[int, int] = None
) -> dict:
    """
    Implementa o algoritmo de Smith-Waterman para alinhamento local de sequências.

    Parâmetros:
        seq_1 (str): Primeira sequência.
        seq_2 (str): Segunda sequência.
        start_pos (tuple[int, int], opcional): Posição inicial (linha, coluna) para iniciar o traceback.
            Se None, utiliza a posição de maior score.

    Retorna:
        dict: Dicionário com as subsequências alinhadas, score, identidade e matches.
    """
    w = Score.GAP.value

    M = len(seq_1)
    N = len(seq_2)
    matrix = np.zeros((M+1, N+1))
    traceback = np.zeros((M+1, N+1), dtype=int)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, M+1):
        for j in range(1, N+1):
            score = get_score(seq_1[i-1], seq_2[j-1])
            diag = matrix[i-1][j-1] + score
            up = matrix[i-1][j] + w
            left = matrix[i][j-1] + w
            vals = [0, up, left, diag]
            best = max(vals)
            matrix[i][j] = best
            # Prioridade: diag > up > left > zero
            if best == diag:
                traceback[i][j] = DIAG
            elif best == up:
                traceback[i][j] = UP
            elif best == left:
                traceback[i][j] = LEFT
            else:
                traceback[i][j] = -1  # zero

            if best > max_score:
                max_score = best
                max_pos = (i, j)

    # Traceback
    aligned_seq_1 = []
    aligned_seq_2 = []
    if start_pos is not None:
        i, j = start_pos
    else:
        i, j = max_pos

    while i > 0 and j > 0 and matrix[i][j] > 0:
        if traceback[i][j] == DIAG:
            aligned_seq_1.append(seq_1[i-1])
            aligned_seq_2.append(seq_2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == UP:
            aligned_seq_1.append(seq_1[i-1])
            aligned_seq_2.append('-')
            i -= 1
        elif traceback[i][j] == LEFT:
            aligned_seq_1.append('-')
            aligned_seq_2.append(seq_2[j-1])
            j -= 1
        else:
            break

    aligned_seq_1 = ''.join(reversed(aligned_seq_1))
    aligned_seq_2 = ''.join(reversed(aligned_seq_2))
    final_score = int(matrix[max_pos[0]][max_pos[1]] if start_pos is None else matrix[start_pos[0]][start_pos[1]])

    # Calcula a porcentagem de identidade
    matches = sum(a == b for a, b in zip(aligned_seq_1, aligned_seq_2))
    identity = matches / len(aligned_seq_1) if aligned_seq_1 else 0

    result = {
        "seq_1": aligned_seq_1,
        "seq_2": aligned_seq_2,
        "score": final_score,
        "identity": identity,
        "match": matches
    }
    return result

def main(id_1: str, id_2: str, start_pos: tuple[int, int] = None):
        seq1 = read_sequence(id_1)
        seq2 = read_sequence(id_2)
        return smith_waterman(seq1, seq2, start_pos)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        res = main(sys.argv[1], sys.argv[2])
        print(res)
    elif len(sys.argv) == 5:
        start_pos = (int(sys.argv[3]), int(sys.argv[4]))
        res = main(sys.argv[1], sys.argv[2], start_pos)
        print(res)
    else:
        print("Uso: python local_alignment.py <tax_id_1> <tax_id_2> [<start_row> <start_col>]")