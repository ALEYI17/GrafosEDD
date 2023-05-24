#include <iostream>
#include <vector>
#include <cstdlib>
#include <stack>
#include <queue>
#include <limits>    // Add this header for std::numeric_limits
#include <algorithm>

template <class T, class U>
class Grafo {
protected:
    std::vector<T> vertices;
    std::vector<std::vector<U>> aristas;
    std::vector<std::vector<U>> aristasTranspuesta;

public:
    Grafo() {}

    void setVertices(const std::vector<T>& v) {
        vertices = v;
    }

    void setAristas(const std::vector<std::vector<U>>& a) {
        aristas = a;
    }

    std::vector<T> getVertices() {
        return vertices;
    }

    std::vector<std::vector<U>> getAristas() {
        return aristas;
    }

    int cantVertices() {
        return vertices.size();
    }

    int cantAristas() {
        int suma = 0;
        for (int i = 0; i < cantVertices(); i++) {
            for (int j = 0; j < cantVertices(); j++) {
                if (aristas[i][j] != 0)
                    suma++;
            }
        }
        return suma;
    }

    int buscarVertice(T vert) {
        for (int i = 0; i < cantVertices(); i++) {
            if (vertices[i] == vert)
                return i;
        }
        return -1;
    }

    bool insertarVertice(T vert) {
        if (buscarVertice(vert) == -1) {
            vertices.push_back(vert);
            aristas.push_back(std::vector<U>(cantVertices(), 0));
            for (int i = 0; i < cantVertices() - 1; i++) {
                aristas[i].push_back(0);
            }
            return true;
        }
        return false;
    }

    bool insertarArista(T ori, T des, U cos) {
        int i_ori = buscarVertice(ori);
        int i_des = buscarVertice(des);
        if (i_ori != -1 && i_des != -1) {
            if (aristas[i_ori][i_des] == 0) {
                aristas[i_ori][i_des] = cos;
                return true;
            }
        }
        return false;
    }

    U buscarArista(T ori, T des) {
        int i_ori = buscarVertice(ori);
        int i_des = buscarVertice(des);
        if (i_ori != -1 && i_des != -1) {
            return aristas[i_ori][i_des];
        }
        return -1;
    }

    bool eliminarVertice(T vert) {
        int indice = buscarVertice(vert);
        if (indice != -1) {
            vertices.erase(vertices.begin() + indice);
            aristas.erase(aristas.begin() + indice);
            for (int i = 0; i < cantVertices(); i++) {
                aristas[i].erase(aristas[i].begin() + indice);
            }
            return true;
        }
        return false;
    }

    bool eliminarArista(T ori, T des) {
        int i_ori = buscarVertice(ori);
        int i_des = buscarVertice(des);
        if (i_ori != -1 && i_des != -1) {
            aristas[i_ori][i_des] = 0;
            return true;
        }
        return false;
    }

    void imprimirMatriz() const {
            for (const auto& row : aristas) {
                for (const auto& val : row) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
    }

    void recorridoPlano(){
        for(const auto& lis : vertices){
            std::cout << lis << " "; 
        }
    }

    void DFS(T inicio) {
        std::vector<bool> visitado(vertices.size(), false); // Track visited vertices
        std::stack<T> pila; // Stack to perform DFS
        int inicioIndex = buscarVertice(inicio);

        if (inicioIndex == -1) {
            std::cout << "El vértice de inicio no existe." << std::endl;
            return;
        }

        pila.push(vertices[inicioIndex]); // Push the starting vertex to the stack
        visitado[inicioIndex] = true; // Mark the starting vertex as visited

        std::cout << "DFS traversal comenzando desde " << inicio << ": ";

        while (!pila.empty()) {
            T actual = pila.top(); // Get the top vertex from the stack
            pila.pop();

            int actualIndex = buscarVertice(actual);

            std::cout << actual << " ";

            // Traverse all adjacent vertices of the current vertex
            for (int i = 0; i < vertices.size(); i++) {
                if (aristas[actualIndex][i] != 0 && !visitado[i]) {
                    pila.push(vertices[i]); // Push the unvisited adjacent vertex to the stack
                    visitado[i] = true; // Mark the adjacent vertex as visited
                }
            }
        }

        std::cout << std::endl;
    }

    void BFS(T inicio){
        std::vector<bool> visitado(vertices.size(), false);
        std::queue<T> cola;
        int inicioIndex = buscarVertice(inicio);

        if (inicioIndex == -1) {
            std::cout << "El vértice de inicio no existe." << std::endl;
            return;
        }
        cola.push(vertices[inicioIndex]);
        visitado[inicioIndex] == true;
        std::cout << "BFS comenzando desde " << inicio << ": ";
        while (!cola.empty()){
            T actual =cola.front();
            cola.pop();
            int actualIndex = buscarVertice(actual);
            std::cout << actual << " ";
            for (int i = 0; i < vertices.size(); i++){
                if (aristas[actualIndex][i] != 0 && !visitado[i]){
                    cola.push(vertices[i]);
                    visitado[i] = true;
                }
            }

        }
        std::cout << std::endl;

    }

     void DFS(int v, std::vector<bool>& visitado, std::stack<T>& pila) {
        visitado[v] = true;

        int vIndex = buscarVertice(vertices[v]);

        for (int i = 0; i < vertices.size(); i++) {
            if (aristas[vIndex][i] != 0 && !visitado[i]) {
                DFS(i, visitado, pila);
            }
        }

        pila.push(vertices[v]);
    }

    void DFS_Transpose(int v, std::vector<bool>& visitado, std::vector<T>& componente) {
        visitado[v] = true;
        componente.push_back(vertices[v]);

        int vIndex = buscarVertice(vertices[v]);

        for (int i = 0; i < vertices.size(); i++) {
            if (aristasTranspuesta[vIndex][i] != 0 && !visitado[i]) {
                DFS_Transpose(i, visitado, componente);
            }
        }
    }

    std::vector<std::vector<T>> stronglyConnectedComponents() {
        std::vector<bool> visitado(vertices.size(), false);
        std::stack<T> pila;

        // Perform DFS on the original graph and fill the stack
        for (int i = 0; i < vertices.size(); i++) {
            if (!visitado[i]) {
                DFS(i, visitado, pila);
            }
        }

        // Transpose the adjacency matrix
        aristasTranspuesta.resize(vertices.size(), std::vector<U>(vertices.size()));
        for (int i = 0; i < vertices.size(); i++) {
            for (int j = 0; j < vertices.size(); j++) {
                aristasTranspuesta[i][j] = aristas[j][i];
            }
        }

        // Reset the visited array
        std::fill(visitado.begin(), visitado.end(), false);

        std::vector<std::vector<T>> componentes;

        // Process the vertices in the stack order
        while (!pila.empty()) {
            T actual = pila.top();
            pila.pop();

            int actualIndex = buscarVertice(actual);

            if (!visitado[actualIndex]) {
                std::vector<T> componente;
                DFS_Transpose(actualIndex, visitado, componente);
                componentes.push_back(componente);
            }
        }

        return componentes;
    }
    std::vector<std::pair<T, T>> getBridgeEdges() {
    std::vector<bool> visited(vertices.size(), false);
    std::vector<int> disc(vertices.size(), -1);
    std::vector<int> low(vertices.size(), -1);
    std::vector<int> parent(vertices.size(), -1);
    std::vector<std::pair<T, T>> bridges;

    int time = 0; // Variable to track discovery time

    // Perform DFS to find bridge edges
    for (int v = 0; v < vertices.size(); v++) {
        if (!visited[v]) {
            getBridgeEdgesUtil(v, visited, disc, low, parent, bridges, time);
        }
    }

    return bridges;
    }

    void getBridgeEdgesUtil(int v, std::vector<bool>& visited, std::vector<int>& disc,
                            std::vector<int>& low, std::vector<int>& parent,
                            std::vector<std::pair<T, T>>& bridges, int& time) {
        visited[v] = true;
        disc[v] = time;
        low[v] = time;
        time++;

        for (int i = 0; i < vertices.size(); i++) {
            if (aristas[v][i] != 0) {
                int u = i;

                if (!visited[u]) {
                    parent[u] = v;
                    getBridgeEdgesUtil(u, visited, disc, low, parent, bridges, time);

                    low[v] = std::min(low[v], low[u]);

                    if (low[u] > disc[v]) {
                        bridges.emplace_back(vertices[v], vertices[u]);
                    }
                } else if (u != parent[v]) {
                    low[v] = std::min(low[v], disc[u]);
                }
            }
        }
    }

    bool isEulerian() const {
    // Check if all vertices have even degrees
    for (int i = 0; i < vertices.size(); i++) {
        int degree = 0;
        for (int j = 0; j < vertices.size(); j++) {
            if (aristas[i][j] != 0)
                degree++;
        }
        if (degree % 2 != 0) {
            return false;
        }
    }
    return true;
    }

    bool isEulerian(const std::vector<std::vector<U>>& aristasCopy) const {
    // Check if all vertices have even degrees
    for (int i = 0; i < vertices.size(); i++) {
        int degree = 0;
        for (int j = 0; j < vertices.size(); j++) {
            if (aristasCopy[i][j] != 0)
                degree++;
        }
        if (degree % 2 != 0) {
            return false;
        }
    }
    return true;
    }


    std::vector<T> fleuryEulerPath() {
        if (!isEulerian()) {
            std::cout << "The graph doesn't have an Eulerian path." << std::endl;
            return std::vector<T>();
        }

        // Make a copy of the adjacency matrix to track the visited edges
        std::vector<std::vector<U>> aristasCopy = aristas;

        // Starting vertex for the Eulerian path
        T startVertex = vertices[0];
        int startVertexIndex = buscarVertice(startVertex);

        std::vector<T> path; // Stores the vertices in the Eulerian path
        path.push_back(startVertex);

        while (true) {
            bool removed = false;

            // Traverse all adjacent vertices of the current vertex
            for (int i = 0; i < vertices.size(); i++) {
                if (aristasCopy[startVertexIndex][i] != 0) {
                    // Remove the edge temporarily
                    U weight = aristasCopy[startVertexIndex][i];
                    aristasCopy[startVertexIndex][i] = aristasCopy[i][startVertexIndex] = 0;

                    // Check if the graph becomes non-Eulerian after removing the edge
                    if (!isEulerian(aristasCopy)) {
                        // Restore the removed edge
                        aristasCopy[startVertexIndex][i] = aristasCopy[i][startVertexIndex] = weight;
                    }
                    else {
                        // Edge can be removed, update the current vertex and add it to the path
                        startVertex = vertices[i];
                        startVertexIndex = buscarVertice(startVertex);
                        path.push_back(startVertex);
                        removed = true;
                        break;
                    }
                }
            }

            // If no edge can be removed, we have the Eulerian path
            if (!removed)
                break;
        }

        return path;
    }

    std::vector<std::pair<T, T>> primMST() {
    std::vector<std::pair<T, T>> mst; // Minimum Spanning Tree
    std::vector<bool> visitado(vertices.size(), false); // Track visited vertices
    std::vector<U> dist(vertices.size(), std::numeric_limits<U>::max()); // Distance values used to pick minimum weight edge
    std::vector<T> parent(vertices.size()); // Array to store constructed MST

    // Set the distance of the first vertex to 0
    dist[0] = 0;
    parent[0] = vertices[0];

    for (int i = 0; i < vertices.size() - 1; i++) {
        // Find the vertex with the minimum distance value
        int minIndex = -1;
        U minDist = std::numeric_limits<U>::max();
        for (int j = 0; j < vertices.size(); j++) {
            if (!visitado[j] && dist[j] < minDist) {
                minDist = dist[j];
                minIndex = j;
            }
        }

        // Mark the picked vertex as visited
        visitado[minIndex] = true;

        // Update the distance values of the adjacent vertices
        for (int j = 0; j < vertices.size(); j++) {
            if (!visitado[j] && aristas[minIndex][j] != 0 && aristas[minIndex][j] < dist[j]) {
                dist[j] = aristas[minIndex][j];
                parent[j] = vertices[minIndex];
            }
        }
    }

    // Construct the minimum spanning tree
    for (int i = 1; i < vertices.size(); i++) {
        mst.push_back(std::make_pair(parent[i], vertices[i]));
    }

    return mst;
    }

    std::pair<std::vector<U>, std::vector<std::vector<int>>> dijkstra(T origen) {
        int inicioIndex = buscarVertice(origen);

        if (inicioIndex == -1) {
            std::cout << "El vértice de origen no existe." << std::endl;
            return std::make_pair(std::vector<U>(), std::vector<std::vector<int>>());
        }

        std::vector<U> distancias(vertices.size(), std::numeric_limits<U>::max()); // Initialize distances to infinity
        std::vector<bool> visitado(vertices.size(), false); // Track visited vertices
        std::vector<std::vector<int>> rutas(vertices.size()); // Track shortest paths

        distancias[inicioIndex] = 0; // Distance from the source vertex to itself is 0
        rutas[inicioIndex].push_back(inicioIndex); // The path to the source vertex only contains itself

        for (int i = 0; i < vertices.size() - 1; i++) {
            int minDistIndex = -1;
            U minDist = std::numeric_limits<U>::max();

            // Find the vertex with the minimum distance
            for (int j = 0; j < vertices.size(); j++) {
                if (!visitado[j] && distancias[j] < minDist) {
                    minDist = distancias[j];
                    minDistIndex = j;
                }
            }

            // Mark the selected vertex as visited
            visitado[minDistIndex] = true;

            // Update the distances and paths of the adjacent vertices
            for (int j = 0; j < vertices.size(); j++) {
                if (!visitado[j] && aristas[minDistIndex][j] != 0) {
                    U newDist = distancias[minDistIndex] + aristas[minDistIndex][j];
                    if (newDist < distancias[j]) {
                        distancias[j] = newDist;
                        rutas[j] = rutas[minDistIndex];
                        rutas[j].push_back(j);
                    }
                }
            }
        }

        return std::make_pair(distancias, rutas);
    }

    std::vector<std::vector<U>> floydWarshall() {
        int numVertices = vertices.size();
        
        std::vector<std::vector<U>> distancias(numVertices, std::vector<U>(numVertices, std::numeric_limits<U>::max()));

        // Initialize the diagonal elements to 0 (distance from a vertex to itself)
        for (int i = 0; i < numVertices; i++) {
            distancias[i][i] = 0;
        }

        // Fill the distances based on the existing edges
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                if (aristas[i][j] != 0) {
                    distancias[i][j] = aristas[i][j];
                }
            }
        }

        // Compute shortest distances using the Floyd-Warshall algorithm
        for (int k = 0; k < numVertices; k++) {
            for (int i = 0; i < numVertices; i++) {
                for (int j = 0; j < numVertices; j++) {
                    if (distancias[i][k] != std::numeric_limits<U>::max() && distancias[k][j] != std::numeric_limits<U>::max()) {
                        U newDist = distancias[i][k] + distancias[k][j];
                        if (newDist < distancias[i][j]) {
                            distancias[i][j] = newDist;
                        }
                    }
                }
            }
        }

        return distancias;
    }
};

int main() {
    // Create a Grafo object
    Grafo<int, int> grafo;

    grafo.insertarVertice(1);
    grafo.insertarVertice(2);
    grafo.insertarVertice(3);
    grafo.insertarVertice(4);
    grafo.insertarVertice(5);


    grafo.insertarArista(1, 4, 3);
    grafo.insertarArista(1, 3, 6);
    grafo.insertarArista(2, 1, 3);
    grafo.insertarArista(3, 4, 2);
    grafo.insertarArista(4, 3, 1);
    grafo.insertarArista(4, 2, 1);
    grafo.insertarArista(5, 2, 4);
    grafo.insertarArista(5, 4, 2);


    grafo.imprimirMatriz();

    std::vector<std::vector<int>> result = grafo.floydWarshall();

    // Print the shortest distances
    for (int i = 0; i < result.size(); i++) {
        for (int j = 0; j < result[i].size(); j++) {
            std::cout << "Shortest distance from vertex " << i << " to vertex " << j << ": " << result[i][j] << std::endl;
        }
    }
    return 0;
}
//Falta
//recorrido en profundidad 
//recorrido en anchura
//componentes conectados (retorne una una multilista con los componentes conectados)
//aristas puente (retorne un vector que diga para cada arista si es o no puente)
//algoritmo de Fleury (retorne la secuencia de indices del camino de Euler)
// algoritmo de prim
// algoritmo de dijkstra