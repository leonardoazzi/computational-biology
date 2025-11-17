# importando as bibliotecas necessárias
import pandas as pd # manipulação de datasets
import numpy as np  # operações matemáticas 

# sklearn é uma biblioteca de data science e machine learning
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

# biblioteca para visualização
import matplotlib.pyplot as plt

# definir o nome do arquivo de dados que usaremos
file_name = 'Prostate_GSE6919_U95B.csv' # https://sbcb.inf.ufrgs.br/cumida

if __name__ == '__main__': 

    # tutorial do pandas: https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html

    # usando a biblioteca pandas para abrir o arquivo .csv com os dados de expressão gênica
    data = pd.read_csv(file_name, delimiter=',', header=0, index_col=0) 
    # exibe um resumo do conjunto de dados
    print(data) 
    # para acessar uma coluna do conjunto de dados, por exemplo a coluna com a condição das amostras (label 'type' nesse arquivo de exemplo):
    print(data['type'])
    # para obtermos uma lista com os identificadores das classes:
    classes = data['type'].unique()
    print(classes)
    # para acessar apenas as amostras saudáveis nesse arquivo de exemplo:
    healthy = data[data['type'] == 'normal']
    print(healthy)
    # para obtermos a média do valor de expressão gênica de cada gene (coluna):
    avg = data.mean()
    print(avg)

    # para separar os dados entre as expressões em X e as classes em Y:
    Y = data['type']
    X = data.drop('type', axis=1)
    print(X)
    print(Y)

    # para melhores resultados pode ser uma boa ideia calcularmos o z-score das expressões antes de continuarmos: https://www.statology.org/z-score-python/

    # lembre-se da necessidade de criar os conjuntos de treinamento, validação e teste: 
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    # https://realpython.com/train-test-split-python-data/

    # ler documentação de SVM: https://scikit-learn.org/stable/modules/svm.html
    clf = svm.SVC() # cria um modelo simples de SVM
    clf.fit(X, Y) # treina o SVM nos pares de entrada X e Y
    # para predizermos a classe de uma amostra desconhecida a, usamos o comando abaixo, onde a é o vetor com a expressão gênica da amostra
    # clf.predict(a)
    prediction = clf.predict(X)
    print(prediction) # mostra as classes preditas para X pelo modelo em clf. Como esse resultado se compara às classes reais em Y?
    # para usar o kernel trick devemos passar o kernel como parâmetro na criação do SVM. Um kernel popular é o rbf (Radial Basis Function)
    # rbf_svc = svm.SVC(kernel='rbf')
    # rbf_svc.fit(X, Y)

    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html?highlight=confusion%20matrix#sklearn.metrics.confusion_matrix
    cm = confusion_matrix(Y, prediction, labels=classes) # cria a matriz de confusão a partir das classes corretas (Y), classes preditas (prediction) e nome das classes
    print(cm)

    # ler documentação de k-means: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html

    kmeans = KMeans(n_clusters=2, random_state=0).fit(X) # cria um modelo simples de k-means com 2 clusters, obtido dos dados em X
    clusters = kmeans.predict(X) 
    print(clusters) # mostra a qual cluster cada amostra foi atribuida, os clusters serão identificados por números inteiros a partir de zero
    # o que acontece se aumentarmos o número de clusters (em n_clusters=2)

    # ler documentação de métricas: https://scikit-learn.org/stable/modules/classes.html?highlight=metrics#module-sklearn.metrics
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html?highlight=silhouette#sklearn.metrics.silhouette_score
    # algumas métricas interessantes: accuracy, sensitivy, specificity, matriz de confusão, Silhouette Coefficient (para clusters)
    # que análises são possíveis a partir dos resultados obtidos?

    # extra: visualizar os dados com PCA
    # https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    # https://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_iris.html
