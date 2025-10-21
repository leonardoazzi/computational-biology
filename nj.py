import sys
import numpy as np
import newick as nw
from ete3 import Tree

# Algoritmo: Neighbour-joining
#
# Entradas:
#   D: matriz de distâncias
#   L: labels indexadas pelas linhas e colunas
#
# Saída:
#   N: Árvore Newick nomeada pelas labels, com distâncias
#
# 0. Inicializar nodos para cada OTU
#   n: número de OTUs ou linhas/colunas em D ou tamanho da lista labels
#
# ENQUANTO n > 2 FAÇA:
#   1. Calcular para cada OTU
#     u_i = sum from (j: j != i) to n of D_ij / (n-2), 
#         onde n é o número de OTUs e i é o índice das OTUs
#   
#   2. Encontrar e escolher i e j com os menores valores de
#     D_ij - u_i - u_j 
#         onde i != j
#   
#   3. Criar um novo nodo (ij) que conecta as OTUs i e j
#     3.1 Calcular o tamanho da aresta entre i e o novo nodo (ij) com:
#         v_i = 1/2 * D_ij + 1/2 * (u_i - u_j)
#     3.2 Calcular o tamanho da aresta entre j e o novo nodo (ij) com:
#         v_j = 1/2 * D_ij + 1/2 * (u_j - u_i)
#   
#   4. Estimar a distância entre o novo nodo (ij) e todo os outros taxons
#     D_(ij)k = (D_ik + D_jk - D_ij) / 2
#         onde k != i != j
#   
#   5. Remover linhas e colunas i e j da matriz D e adicionar linha e coluna para (ij)
#     5.1 Remover labels de i e j, e adicionar label de ij
#
#   6. Cria aresta entre ( ij ) e os nodos i e j, com distâncias v_i e v_j
#
# FIM
# 
# 7. Criar nodo (lm) conectando os últimos nodos, e atribuindo distância D_lm
#
# Retornar árvore Newick
# 
# FIM

class Newick:
    def __init__(self, matrix: np.ndarray, labels: list):
        """Cria um nodo nomeado para cada label.

        Args:
            matrix (np.matrix): matriz de distâncias
            labels (list): nomes das OTUs, indexado pelas linhas e colunas de d_matrix
        """        
        if matrix.shape[0] != len(labels) or matrix.shape[0] != len(labels):
            raise Exception("O número de labels é diferente do número de linhas/colunas.")
        
        self.nodes = [nw.Node.create() for _ in labels]
        for idx, node in enumerate(self.nodes):
            node.name = f"{labels[idx]}"
        
    def create_node(self, name: str):
        """Retorna um novo nodo nomeado

        Args:
            name (str): nome do nodo

        Returns:
            nw.Node: objeto da classe nw.Node
        """        
        new_node = nw.Node.create()
        new_node.name = f"{name}"

        return new_node

def get_matrix(matrix_file: str, labels_file: str) -> tuple[np.ndarray, list]:
    """Retorna uma tupla com a matriz de distâncias a partir de um arquivo .txt e as labels a partir de um .txt separado por vírgulas.

    Args:
        matrix_file (str): caminho para o arquivo da matriz em .txt
        labels_file (str): caminho para o arquivo das labels em .txt

    Returns:
        tuple[np.ndarray, list]: tupla com (matrix, labels)
    """    
    with open(labels_file, 'r') as f:
        labels_str = f.read()
    
    labels = [x.strip() for x in labels_str.split(',')]

    matrix = np.loadtxt(matrix_file)

    return matrix, labels

def compute_u(matrix: np.ndarray, clusters: list):
    """Calcula para cada OTU
        u_i = sum from (j: j != i) to n of D_ij / (n-2)
    onde n é o número de OTUs e i é o índice das OTUs

    Args:
        matrix (np.ndarray): matriz de distâncias D
        tree (Newick): objeto da classe Newick

    Returns:
        list: lista com valores de u_i, onde i representa cada OTU
    """    
    n = len(clusters)
    u_list = []
    
    for i, _ in enumerate(matrix[0]):
        # Somatório
        acc = 0
        for j, _ in enumerate(matrix[1]):
            if j != i:
                acc = acc + matrix[i][j] / (n - 2)
        u_list.append(acc)
    
    return u_list

def find_min_ij(matrix: np.ndarray, u_list: list)-> tuple[int, int]:
    """Encontrar e escolher i e j com os menores valores de
        D_ij - u_i - u_j 
    onde i != j

    Args:
        matrix (np.ndarray): matriz de distâncias D
        u_list (list): lista de valores obtidos em compute_u()

    Returns:
        tuple[int, int]: tupla com índices (i, j)
    """    
#   2. 

    u = np.array(u_list, dtype=float)

    # Utiliza broadcasting para computar D_ij - u_i - u_j diretamente
    D_minus_u = matrix - u[:, None] - u[None, :]

    np.fill_diagonal(D_minus_u, np.inf) # Atribui infinito à diagonal i=j
    
    # Encontra índice com menores valores
    i,j = np.unravel_index(np.argmin(D_minus_u), D_minus_u.shape) 
    
    return (i, j)

def compute_edges(matrix: np.ndarray, u_list: list, i: int, j: int) -> tuple[float, float]:
    """_summary_

    Args:
        matrix (np.ndarray): matriz de distâncias D
        u_list (list): lista de valores obtidos em compute_u()
        i (int): índice para o menor valor obtido em find_min_ij()
        j (int): índice para o menor valor obtido em find_min_ij()

    Returns:
        tuple[float, float]: _description_
    """
    v_i = 1/2 * matrix[i,j] + 1/2 * (u_list[i] - u_list[j])
    v_j = 1/2 * matrix[i,j] + 1/2 * (u_list[j] - u_list[i])

    return v_i, v_j

def compute_distances(matrix: np.ndarray, i: int, j: int) -> np.ndarray:
    """Estimar a distância entre o novo nodo (ij) e todo os outros taxons
        D_(ij)k = (D_ik + D_jk - D_ij) / 2
    onde k != i != j

    Args:
        matrix (np.ndarray): matriz de distâncias D
        i (int): índice para o menor valor obtido em find_min_ij()
        j (int): índice para o menor valor obtido em find_min_ij()

    Returns:
        np.ndarray: vetor D_ijk
    """    
    # Cria uma máscara para obter os índices diferentes de i e j
    n = matrix.shape[0]
    mask = np.ones(n, dtype=bool)
    mask[[i,j]] = False
    k_indices = np.arange(n)[mask]

    # Computa D_ijk com broadcasting, gerando um array com elementos ij_k
    D_ijk = (matrix[i, k_indices] + matrix[j, k_indices] - matrix[i, j]) / 2

    return D_ijk

def update_distances(matrix: np.ndarray, D_ijk: np.ndarray, i:int, j:int) -> np.ndarray:
    """Remover linhas e colunas i e j da matriz D e adicionar linha e coluna para D_(ij)k

    Args:
        matrix (np.ndarray): matriz de distâncias D
        D_ijk (np.ndarray): vetor D_ijk
        i (int): índice para o menor valor obtido em find_min_ij()
        j (int): índice para o menor valor obtido em find_min_ij()

    Returns:
        np.ndarray: matriz atualizada
    """    
    updated_matrix = matrix.copy()

    updated_matrix = np.delete(updated_matrix, [i, j], axis=0)
    updated_matrix = np.delete(updated_matrix, [i, j], axis=1)

    n = updated_matrix.shape[0]
    new_matrix = np.zeros((n+1, n+1))

    # Adiciona zero da diagonal
    D_ijk = np.append(D_ijk, 0.0)

    # Copia a matriz original até as linhas e colunas n
    new_matrix[:n, :n] = updated_matrix

    new_matrix[-1, :n+1] = D_ijk # última linha
    new_matrix[:n+1, -1] = D_ijk # última coluna
    
    return new_matrix

def run(d_matrix: np.ndarray, labels: list):
    """Retorna a árvore filogenética obtida com o algoritmo Neighbour-joining.

    Args:
        d_matrix (np.ndarray): matriz de distâncias
        labels (list): nomes das OTUs, indexado pelas linhas e colunas de d_matrix

    Returns:
        tuple: árvore filogenética em formato Newick (sem distâncias)
    """    
    
    # 0. Inicializar nodos para cada OTU
    #   n: número de OTUs ou linhas/colunas em D ou tamanho da lista labels
    #       n = len(Newick.nodes)
    tree = Newick(d_matrix, labels)

    clusters = tree.nodes.copy()
    
    while len(clusters) > 2:
    
        # 1. Calcular u_i para cada OTU
        u_list = compute_u(d_matrix, clusters)
        print("1. u_list:", u_list)

        # 2. Encontrar i e j com os menores valores D_ij - u_i - u_j
        i, j = find_min_ij(d_matrix, u_list)
        print("2. i,j:", i, j)
        
        # 3. Criar novo nodo (ij) que conecta OTUs i e j
        i_name = clusters[i].name
        j_name = clusters[j].name
        node_ij = tree.create_node("")
        print("3. Novo nodo:", node_ij)

        # 3.1. Calcular o tamanho da aresta entre os nodos i e j com o nodo (ij)
        v_i, v_j = compute_edges(d_matrix, u_list, i, j)
        
        # 4. Estimar distâncias entre nodo (ij) e os k restantes
        D_ijk = compute_distances(d_matrix, i, j)
        print("4. D_ijk", D_ijk)

        # 5. Atualizar matriz D removendo i e j, e adicionando D_ijk
        d_matrix = update_distances(tree, d_matrix, D_ijk, i, j)

        # 5.1 Atualizando índice de nodos
        node_i = clusters[i]
        node_j = clusters[j]

        ## Remove nodos i e j dos clusters
        clusters = [v for idx, v in enumerate(clusters) if idx not in [i,j]]

        # 6. Cria aresta entre ( ij ) e os nodos i e j
        node_ij.add_descendant(node_i)
        node_ij.add_descendant(node_j)

        # Atualiza distâncias de i e de j para (ij)
        node_i.length = v_i
        node_j.length = v_j

        ## Adiciona nodo ij aos clusters e à árvore
        clusters.append(node_ij)
        tree.nodes.append(node_ij)

        print(i,j)
        print("clusters", clusters)
        print("nodes", tree.nodes)
        # print(nw.dumps(node_ij))
        print(node_ij.ascii_art())
        print(d_matrix)
    
    # 7. Criar nodo lm conectando os últimos nodos l e m, e atribuindo distância D_lm
    print("clusters", clusters)
    node_l = clusters[0]
    node_m = clusters[1]

    node_m.add_descendant(node_l)
    print(node_m.ascii_art())
    print(nw.dumps(node_ij))
    newick_tree = nw.dumps(node_ij)
    t = Tree(newick_tree, format=0)

    t.show()

if __name__ == "__main__":
    if len(sys.argv) == 3:
        matrix_file = sys.argv[1]
        labels_file = sys.argv[2]
    else:
        print("Uso: python nj.py <matrix_file> <labels_file>")

    matrix, labels = get_matrix(matrix_file, labels_file)
    newick = run(matrix, labels)
    # print(newick)