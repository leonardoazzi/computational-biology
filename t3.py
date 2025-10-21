from newick import Node, dumps

teste = Node.create()
teste.name = "A"
teste2 = Node.create()
teste2.name = "B" 


teste3 = Node.create()
teste3.name = "C"

teste3.add_descendant(teste2)
teste3.add_descendant(teste)
teste.length = 3
teste2.length = 2

print(teste3.ascii_art())

print(dumps(teste3))