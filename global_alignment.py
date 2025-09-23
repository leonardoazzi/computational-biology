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

def needleman_wunsch(seq_1: str, seq_2: str) -> tuple[str, str, int, float]:
    """
    Implementa o algoritmo de Needleman-Wunsch para alinhamento global de sequências.

    Parâmetros:
        seq_1 (str): Primeira sequência a ser alinhada.
        seq_2 (str): Segunda sequência a ser alinhada.

    Retorna:
        tuple[str, str, int, float]: Tupla contendo as versões alinhadas de seq_1 e seq_2,
        a pontuação final e a porcentagem de identidade.
    """
    w = Score.GAP.value

    M = len(seq_1)
    N = len(seq_2)
    matrix = np.zeros((M+1, N+1))

    # Inicialização da primeira coluna (percorrendo linhas)
    for i in range(M+1):
        matrix[i][0] = w * i

    # Inicialização da primeira linha (percorrendo colunas)
    for j in range(N+1):
        matrix[0][j] = w * j

    for i in range(1, M+1):
        for j in range(1, N+1):
            score = get_score(seq_1[i-1], seq_2[j-1])
            diag = matrix[i-1][j-1] + score
            up = matrix[i-1][j] + w
            left = matrix[i][j-1] + w
            matrix[i][j] = max(diag, up, left)

    # Traceback
    aligned_seq_1 = []
    aligned_seq_2 = []
    i, j = M, N

    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + get_score(seq_1[i-1], seq_2[j-1]):
            aligned_seq_1.append(seq_1[i-1])
            aligned_seq_2.append(seq_2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + w:
            aligned_seq_1.append(seq_1[i-1])
            aligned_seq_2.append('-')
            i -= 1
        else:
            aligned_seq_1.append('-')
            aligned_seq_2.append(seq_2[j-1])
            j -= 1

    aligned_seq_1 = ''.join(reversed(aligned_seq_1))
    aligned_seq_2 = ''.join(reversed(aligned_seq_2))
    final_score = int(matrix[M][N])

    # Calcula a porcentagem de identidade
    matches = sum(a == b for a, b in zip(aligned_seq_1, aligned_seq_2))
    identity = matches / len(aligned_seq_1)

    # print("Alinhamento global:")
    
    # print(aligned_seq_1)
    # print(aligned_seq_2)
    # print(f"Pontuação final: {final_score}")
    # print(f"Matches: {matches}")
    # print(f"Tamanho: {len(aligned_seq_1)}")
    # print(f"Identidade: {identity:.2f}%\n")

    result = {
        "seq_1": aligned_seq_1,
        "seq_2": aligned_seq_2,
        "score": final_score,
        "identity": identity,
        "match": matches
    }
    return result

def main(id_1: str, id_2: str):
        seq1 = read_sequence(id_1)
        seq2 = read_sequence(id_2)
        return needleman_wunsch(seq1, seq2)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Uso: python global_alignment.py <tax_id_1> <tax_id_2>")