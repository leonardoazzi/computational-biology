import numpy as np
from enum import Enum
import sys
import json
import re

class Score:
    """
    Classe Score para representar os parâmetros de pontuação no alinhamento de sequências.

    Atributos:
        MATCH (int): Pontuação para um match entre dois elementos.
        MISMATCH (int): Penalidade para um mismatch entre dois elementos.
        GAP (int): Penalidade para a introdução de um gap.

    Métodos:
        __init__(match=0, mismatch=0, gap=0):
            Inicializa um objeto Score com valores de match, mismatch e gap especificados.
        from_dict(scoring_dict):
            Método de classe para criar um objeto Score a partir de um dicionário com as chaves 'MATCH', 'MISMATCH' e 'GAP'.
    """
    def __init__(self, match=0, mismatch=0, gap=0):
        self.MATCH = match
        self.MISMATCH = mismatch
        self.GAP = gap

    @classmethod
    def from_dict(cls, scoring_dict: dict) -> 'Score':
        """Cria um objeto Score a partir de um dicionário de pontuação.

        Args:
            scoring_dict (dict): Dicionário contendo as chaves 'MATCH', 'MISMATCH' e 'GAP'.

        Returns:
            Score: Um objeto Score com os valores correspondentes do dicionário.
        """
        return cls(
            match=scoring_dict.get("MATCH"),
            mismatch=scoring_dict.get("MISMATCH"),
            gap=scoring_dict.get("GAP")
        )

def read_sequence(directory_path: str, file_format: str, tax_id: str) -> str:
    with open(f"{directory_path}/dataset/{tax_id}.{file_format}", "r") as file:
        # @TODO: processar arquivos FASTA e arquivos TXT
        sequence = re.sub(r'^>.*\n', '', file.read().strip())
        sequences = []
        sequences = [seq.replace('\n', '').strip() for seq in sequences if seq.strip()]

        return sequence

def read_json(file_path: str) -> dict:
    """Lê um arquivo JSON e o converte em um dicionário.

    Args:
        file_path (str): Caminho para o arquivo JSON.

    Raises:
        ValueError: Se o arquivo não for um JSON válido.

    Returns:
        json_dict: Dicionário contendo os dados do JSON.
    """
    try:
        with open(file_path, "r") as file:
            json_dict = json.load(file)
    except json.JSONDecodeError:
        raise ValueError(f"File '{file_path}' is not a valid JSON file.")
    return json_dict

def get_score(char_1: str, char_2: str, scoring: Score) -> int:
    """Calcula a pontuação entre dois caracteres com base nas regras de pontuação fornecidas.

    Args:
        char_1 (str): primeiro caractere a ser comparado.
        char_2 (str): segundo caractere a ser comparado.
        scoring (Score): objeto de pontuação contendo as regras de match, mismatch e gap.

    Returns:
        int: pontuação calculada entre os dois caracteres.
    """
    score = 0
    if char_1 == char_2:
        score = scoring.MATCH
    else:
        score = scoring.MISMATCH
    
    return score

def needleman_wunsch(seq_1: str, seq_2: str, scoring: Score) -> tuple[str, str, int, float]:
    """
    Implementa o algoritmo de Needleman-Wunsch para alinhamento global de sequências.

    Args:
        seq_1 (str): Primeira sequência a ser alinhada.
        seq_2 (str): Segunda sequência a ser alinhada.
        scoring_dict (dict): Dicionário contendo os valores de pontuação para match, mismatch e gap.

    Returns:
        tuple[str, str, int, float]: Tupla contendo as versões alinhadas de seq_1 e seq_2,
        a pontuação final e a porcentagem de identidade.
    """

    M = len(seq_1)
    N = len(seq_2)
    matrix = np.zeros((M+1, N+1))

    w = scoring.GAP

    # Inicialização da primeira coluna (percorrendo linhas)
    for i in range(M+1):
        matrix[i][0] = w * i

    # Inicialização da primeira linha (percorrendo colunas)
    for j in range(N+1):
        matrix[0][j] = w * j

    for i in range(1, M+1):
        for j in range(1, N+1):
            score = get_score(seq_1[i-1], seq_2[j-1], scoring)
            diag = matrix[i-1][j-1] + score
            up = matrix[i-1][j] + w
            left = matrix[i][j-1] + w
            matrix[i][j] = max(diag, up, left)

    # Traceback
    aligned_seq_1 = []
    aligned_seq_2 = []
    i, j = M, N

    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + get_score(seq_1[i-1], seq_2[j-1], scoring):
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

    result = {
        "seq_1": aligned_seq_1,
        "seq_2": aligned_seq_2,
        "score": final_score,
        "identity": identity,
        "match": matches
    }

    return result

def main(directory_path:str, file_format: str, tax_id_1: str, tax_id_2: str, scoring_filename: str):
    """Realiza o alinhamento global entre duas sequências identificadas por seus tax_ids.

    Args:
        tax_id_1 (str): ID do primeiro organismo.
        tax_id_2 (str): ID do segundo organismo.
        scoring_filepath (str): Caminho para o arquivo JSON de pontuação.

    Returns:
        dict: Resultados do alinhamento global.
    """
    seq1 = read_sequence(directory_path, file_format, tax_id_1)
    seq2 = read_sequence(directory_path, file_format, tax_id_2)

    scoring_filepath = f"{directory_path}/{scoring_filename}.json"
    scoring_dict = read_json(scoring_filepath)
    scoring = Score.from_dict(scoring_dict)

    return needleman_wunsch(seq1, seq2, scoring)

if __name__ == "__main__":
    if len(sys.argv) == 6:
        res = main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        print(res)
    else:
        print("Uso: python global_alignment.py <directory_path> <file_format> <tax_id_1> <tax_id_2> <scoring_filename>")