import sys
import numpy as np

class Cluster:
    def __init__(self, dist_matrix: np.matrix):
        """Cria um cluster para cada linha ou coluna da matriz de distâncias

        Args:
            dist_matrix (np.matrix): matriz de distâncias
        """        
        N = dist_matrix.shape[0]
        self.clusters = [(x) for x in list(range(N))]        

    def group_clusters(self, c_a, c_b):
        """Agrupa dois clusters Newick.

        Args:
            c_a (tuple): primeiro cluster
            c_b (tuple): segundo cluster
        """        
        self.clusters.insert(0, (c_a, c_b)) # Adiciona no início

        a_idx = self.clusters.index(c_a)
        b_idx = self.clusters.index(c_b)

        # Remover ambos os clusters sem alterar os índices antes da remoção
        for idx in sorted([a_idx, b_idx], reverse=True):
            self.clusters.pop(idx)

    def complement_clusters(self, c_a, c_b):
        """Retorna uma lista com os clusters complementares a c_a e c_b.

        Args:
            c_a (tuple): primeiro cluster
            c_b (tuple): segundo cluster

        Returns:
            list: lista com clusters existentes diferentes de c_a e c_b.
        """        
        return [i for _, i in enumerate(self.clusters) if i not in [c_a, c_b]]
    
def run(d_matrix: np.ndarray):
    """Retorna a árvore filogenética obtida com o algoritmo UPGMA.

    Args:
        d_matrix (np.ndarray): matriz de distâncias

    Returns:
        tuple: árvore filogenética em formato Newick (sem distâncias)
    """    
    # 1. Inicialização
    # - Atribuir cada x_i ao seu cluster C_i
    C = Cluster(d_matrix)

    while len(C.clusters) > 1:
        old_cluster = C.clusters.copy()

        # 2. Iteração
        # - Encontrar dois clusters C_i e C_j cujo d_ij é mínimo
        d_matrix_tri = d_matrix.copy()

        # Aplica matriz triangular, com diagonal e triângulo superior em np.inf,
        # para que seja encontrado o valor mínimo do triângulo inferior.
        d_matrix_tri[np.triu_indices_from(d_matrix_tri, k=0)] = np.inf
        d_min = np.unravel_index(np.argmin(d_matrix_tri, axis=None), d_matrix.shape)

        C_i, C_j = d_min # Índices do valor mínimo. Vai indicar quais clusters manipular
        C_comp = C.complement_clusters(C.clusters[C_i], C.clusters[C_j])

        # - C_k = C_i U C_j
        C.group_clusters(C.clusters[C_i], C.clusters[C_j])

        # Criar cópia da matriz original, deletando C_i e C_j
        d_matrix_new = d_matrix.copy()
        d_matrix_new = np.delete(d_matrix_new, C_i, axis=0) # Deletar linhas C_i e C_j
        d_matrix_new = np.delete(d_matrix_new, C_j, axis=0) # Deletar colunas C_i e C_j
        d_matrix_new = np.delete(d_matrix_new, C_i, axis=1) # Deletar linhas C_i e C_j
        d_matrix_new = np.delete(d_matrix_new, C_j, axis=1) # Deletar colunas C_i e C_j

        # Adicionar uma coluna e linha C_k
        d_matrix_new = np.c_[np.ones(d_matrix_new.shape[0]) * np.inf, d_matrix_new]
        d_matrix_new = np.r_[[np.ones(d_matrix_new.shape[1]) * np.inf], d_matrix_new]
        d_matrix_new[0][0] = 0 # mantém diagonal zerada

        # Computar novas distâncias para o cluster C_k
        # - Atualizar matriz com d_kz = (d_iz + d_jz) / 2
        for c_alvo in C_comp:
            # Busca distâncias na matriz original entre C_i e C_j com outros clusters c_k
            c_k = old_cluster.index(c_alvo) # Obtém índice do cluster na matriz original

            d_ik = d_matrix[C_i][c_k]
            d_jk = d_matrix[C_j][c_k]

            d_mk = (d_ik + d_jk) / 2.0

            # Atualiza distâncias na nova matriz
            c_k = C.clusters.index(c_alvo) # Obtém índice do cluster na matriz nova

            d_matrix_new[0][c_k] = d_mk
            d_matrix_new[c_k][0] = d_mk

        # Atualiza a matriz principal
        d_matrix = d_matrix_new

    # 3. Terminação
    # - Quando restar apenas um cluster C_i
    return C.clusters[0]
    
if __name__ == "__main__":
    dir = "./2_phylogenetics"
    matrix = np.loadtxt(f'{dir}/dist_matrix.txt')
    newick = run(matrix)
    print(newick)