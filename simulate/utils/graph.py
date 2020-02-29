class DirectedGraph(object):

    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[] for _ in range(self.V)]

    def add_edge(self, src, dest):
        if dest not in self.graph[src]:
            self.graph[src].append(dest)

    def __repr__(self):
        return_string = ""
        for i in range(len(self.graph)):
            return_string += "{} -> {}\n".format(i, str(self.graph[i]))
        return return_string


class UndirectedGraph(DirectedGraph):

    def __init__(self, vertices):
        super(UndirectedGraph, self).__init__(vertices)

    def add_edge(self, src, dest):
        super(UndirectedGraph, self).add_edge(src, dest)
        super(UndirectedGraph, self).add_edge(dest, src)
