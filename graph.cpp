#include <bits/stdc++.h>
using namespace std;

void printVector(vector<int> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}
// This class is the base class for other graph-related classes such as DirectedGraph, WeightedGraph, and WeightedDirectedGraph
class Graph
{
private:
    // making this private so only friend class i.e. DirectedGraph can access it
    vector<vector<int>> adjList;
    void numConnectedComponentsHelper(int node, vector<int> &visited);
    bool isCyclicHelper(int node, vector<int> &visited, int parent);

protected:
    //  making this protected as numNodes will be common for all graph types
    int numNodes;

public:
    Graph(int nodes); // constructor
    virtual ~Graph() {}
    virtual void addEdge(int v, int w);
    virtual void removeEdge(int v, int w);
    // this display function will be used both in Graph and DirectedGraph classes
    virtual void display();
    virtual vector<int> bfs(int node);
    virtual vector<int> dfs(int node);

    // FOR YOU----> // used call by reference, for just using the address, instead of unnecessary copying
    virtual void dfsHelper(int node, vector<int> &visited, vector<int> &ans);

    bool isCyclic();
    int numConnectedComponents();
    virtual int ShortestPath(int start, int destination);
    // this way only DerictedGraph class can access adjList
    friend class DirectedGraph;
    friend ostream & operator << (ostream &out, const Graph &graph);
};

ostream & operator << (ostream &out, const Graph &graph)
{
    for (int node = 0; node < graph.numNodes; node++){
        out << node << " -> ";
        for(auto neighbour : graph.adjList[node]){
            out << neighbour << " ";
        }
        out << endl;
    }
    return out;
}
// Graph::Graph(int nodes): This is the constructor for the Graph class, which takes an integer argument vertices to specify the number of nodes in the graph.
// : numNodes(nodes): This is a member initialization list, which initializes the numNodes attribute with the value of the vertices argument.
Graph::Graph(int nodes) : numNodes(nodes)
{
    // Resizes the adjacency list to accommodate the specified number of vertices
    adjList.resize(nodes);
}
void Graph::addEdge(int v, int w)
{
    adjList[v].push_back(w);
    adjList[w].push_back(v);
}
void Graph::removeEdge(int v, int w)
{
    adjList[v].erase(find(adjList[v].begin(), adjList[v].end(), w));
    adjList[w].erase(find(adjList[w].begin(), adjList[w].end(), v));
}
void Graph::display()
{
    for (int i = 0; i < numNodes; i++)
    {
        cout << i << " -> ";
        for (int j = 0; j < adjList[i].size(); j++)
        {
            cout << adjList[i][j] << " ";
        }
        cout << endl;
    }
}
vector<int> Graph::bfs(int vertex)
{
    vector<int> ans;
    queue<int> q;
    q.push(vertex);
    vector<int> visited(numNodes, 0);
    visited[vertex] = 1;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : adjList[node])
        {
            if (!visited[neighbour])
            {
                visited[neighbour] = 1;
                q.push(neighbour);
            }
        }
    }
    return ans;
}
void Graph::dfsHelper(int node, vector<int> &visited, vector<int> &ans)
{
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour])
        {
            dfsHelper(neighbour, visited, ans);
        }
    }
}
vector<int> Graph::dfs(int node)
{
    vector<int> ans;
    vector<int> visited(numNodes, 0);
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        dfsHelper(neighbour, visited, ans);
    }
    return ans;
}
bool Graph::isCyclicHelper(int node, vector<int>& visited, int parent)
{
    visited[node] = 1;
    for(auto neighbour : adjList[node]){
        if(!visited[neighbour]){
            if(isCyclicHelper(neighbour, visited, node)){
                return true;
            }
        }else if(neighbour != parent){
            return true;
        }
    }
    return false;
}
bool Graph::isCyclic()
{
    vector<int> visited(numNodes, 0);
    for(int node = 0; node < numNodes; node++){
        if(!visited[node])
        {
            if(isCyclicHelper(node, visited, -1)==true)
            {
                return true;
            }
        }
    }
    return false;
}
void Graph::numConnectedComponentsHelper(int node, vector<int> &visited)
{
    visited[node] = 1;
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour])
        {
            numConnectedComponentsHelper(neighbour, visited);
        }
    }
}
int Graph::numConnectedComponents()
{
    vector<int> visited(numNodes, 0);
    int result = 0;
    for (int node = 0; node < numNodes; node++)
    {
        if (!visited[node])
        {
            numConnectedComponentsHelper(node, visited);
            result++;
        }
    }
    return result;
}
int Graph::ShortestPath(int start, int destination)
{
    // FOR YOU---> // this array is used to store the minimum distance from start to any node encountered during the loop
    vector<int> dp(numNodes, 1e9);
    dp[start] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, start));
    while(!q.empty())
    {
        // lambda capture 
        int distance = q.top().first;
        int node = q.top().second;
        q.pop();
        if (node == destination) break;
        for (auto neighbour : adjList[node])
        {
            if (dp[neighbour] > distance + 1)
            {
                dp[neighbour] = distance + 1;
                q.push(make_pair(dp[neighbour], neighbour));
            }
        }
    }
    // -1 for not reachable from start to destination
    if (dp[destination] == 1e9) return -1;
    return dp[destination];
}
class DirectedGraph : public Graph
{
private:
    bool isCyclicHelper(int node, vector<int>& visited);
public:
    DirectedGraph(int nodes);
    void addEdge(int v, int w);
    void removeEdge(int v, int w);
    // display function of Graph class will be used
    // topological sort using Kahn algorithm
    virtual vector<int> topoSort();
    bool isCyclic();
};
DirectedGraph::DirectedGraph(int nodes) : Graph(nodes) {}
void DirectedGraph::addEdge(int v, int w)
{
    adjList[v].push_back(w); // v:source  w:destination
}
void DirectedGraph::removeEdge(int v, int w)
{
    adjList[v].erase(find(adjList[v].begin(), adjList[v].end(), w));
}

vector<int> DirectedGraph::topoSort()
{
    vector<int> ans;
    vector<int> indegree(numNodes, 0);
    queue<int> q;
    for (int node = 0; node < numNodes; node++)
    {
        for (auto neighbour : adjList[node])
        {
            indegree[neighbour]++;
        }
    }
    for (int node = 0; node < numNodes; node++)
    {
        if (indegree[node] == 0)
        {
            q.push(node);
        }
    }
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : adjList[node])
        {
            indegree[neighbour]--;
            if (indegree[neighbour] == 0)
            {
                q.push(neighbour);
            }
        }
    }
    return ans;
}
bool DirectedGraph::isCyclicHelper(int node, vector<int>& visited)
{
    visited[node] = 1;
    for (auto neighbour : adjList[node]){
        if(!visited[neighbour]){
            if(isCyclicHelper(neighbour, visited)){
                return true;
            }
        }else if (visited[neighbour] == 1){
            return true;
        }
    }
    visited[node] = 2;
    return false;
}
bool DirectedGraph::isCyclic()
{
    vector<int> visited(numNodes, 0);
    for(int node = 0; node < numNodes; node++){
        if (!visited[node])
        {
            if(isCyclicHelper(node, visited)==true)
            {
                return true;
            }
        }
    }
    return false;
}

class WeightedGraph : public Graph
{
    // protected so that WeightedDirectedGraph class can access it
protected:
    vector<vector<pair<int, int>>> adjList;
    void numConnectedComponenetsHelper(int node, vector<int> &visited);
    bool isCyclicHelper(int node, vector<int> & visited, int parent);

public:
    WeightedGraph(int nodes);
    // since number of arguments are different, it is the case of method overloading(consider it different function all together)
    // so we need to again define it to be virtual here(for appropriate binding behaviour in this class and WeightedDirectedGraph class)
    virtual void addEdge(int v, int w, int weight);
    virtual void removeEdge(int v, int w);
    vector<int> bfs(int node);
    vector<int> dfs(int node);
    void dfsHelper(int node, vector<int> &visited, vector<int> &ans);

    int numConnectedComponents();
    virtual int ShortestPath(int start, int target);
    vector<pair<int,int>> MinimumSpanningTree();

    // display has same number of arguments so no need to define it to be virtual here(as already done in Graph class)
    void display();
    friend ostream & operator << (ostream &out, const WeightedGraph &graph);
};

vector<int> WeightedGraph::bfs(int vertex)
{
    vector<int> ans;
    queue<int> q;
    q.push(vertex);
    vector<int> visited(numNodes, 0);
    visited[vertex] = 1;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : adjList[node])
        {
            if (!visited[neighbour.first])
            {
                visited[neighbour.first] = 1;
                q.push(neighbour.first);
            }
        }
    }
    return ans;
}

void WeightedGraph::dfsHelper(int node, vector<int> &visited, vector<int> &ans)
{
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour.first])
        {
            dfsHelper(neighbour.first, visited, ans);
        }
    }
}

vector<int> WeightedGraph::dfs(int node)
{
    vector<int> ans;
    vector<int> visited(numNodes, 0);
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        dfsHelper(neighbour.first, visited, ans);
    }
    return ans;
}

ostream & operator << (ostream &out, const WeightedGraph &graph){
    for (int node = 0; node < graph.numNodes; node++){
        out << node << " -> ";
        for(auto neighbour : graph.adjList[node]){
            out << neighbour.first << " " << neighbour.second << "  ";
        }
        out << endl;
    }
    return out;
}

WeightedGraph::WeightedGraph(int nodes) : Graph(nodes)
{
    // since adjList of this class is different from Graph class, we need to resize it here
    adjList.resize(numNodes);
}
void WeightedGraph::addEdge(int v, int w, int weight)
{
    adjList[v].push_back(make_pair(w, weight));
    adjList[w].push_back(make_pair(v, weight));
}
void WeightedGraph::removeEdge(int v, int w)
{
    for (auto itr = adjList[v].begin(); itr != adjList[v].end(); itr++)
    {
        if (itr->first == w)
        {
            adjList[v].erase(itr);
            break;
        }
    }
    for (auto itr = adjList[w].begin(); itr != adjList[w].end(); itr++)
    {
        if (itr->first == v)
        {
            adjList[w].erase(itr);
            break;
        }
    }
}
void WeightedGraph::numConnectedComponenetsHelper(int node, vector<int> & visited)
{
    visited[node] = 1;
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour.first])
        {
            WeightedGraph::numConnectedComponenetsHelper(neighbour.first, visited);
        }
    }
}
int WeightedGraph::numConnectedComponents()
{
    vector<int> visited(numNodes, 0);
    int result = 0;
    for (int node = 0; node < numNodes; node++)
    {
        if (!visited[node])
        {
            WeightedGraph::numConnectedComponenetsHelper(node, visited);
            result++;
        }
    }
    return result;
}
bool WeightedGraph::isCyclicHelper(int node, vector<int> & visited, int parent)
{
    visited[node] = 1;
    for(auto neighbour : adjList[node]){
        if(!visited[neighbour.first]){
            if(WeightedGraph::isCyclicHelper(neighbour.first, visited, node)){
                return true;
            }
        }else if(neighbour.first != parent){
            return true;
        }
    }
    return false;
}
int WeightedGraph::ShortestPath(int start, int destination)
{
    vector<int> dp(numNodes, 1e9);
    dp[start] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, start));
    while(!q.empty())
    {
        int distance = q.top().first;
        int node = q.top().second;
        q.pop();
        if (node == destination) break;
        for (auto neighbour : adjList[node])
        {
            if (dp[neighbour.first] > distance + neighbour.second)
            {
                dp[neighbour.first] = distance + neighbour.second;
                q.push(make_pair(dp[neighbour.first], neighbour.first));
            }
        }
    }
    // return -1 for not reachable from start to destination
    if (dp[destination] == 1e9) return -1;
    return dp[destination];
}
vector<pair<int,int>> WeightedGraph::MinimumSpanningTree()
{
    vector<int> included(numNodes, 0);
    // stores the minimum edge weight through which the any node was encountered
    vector<int> dp(numNodes, 1e9);
    dp[0] = 0;
    vector<int> parent(numNodes);
    parent[0] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, 0));
    while(!q.empty())
    {
        int distance = q.top().first;
        int node = q.top().second;
        q.pop();
        if (included[node] == 1) continue;
        included[node] = 1;
        for (auto neighbour : adjList[node])
        {
            if (!included[neighbour.first] && dp[neighbour.first] > neighbour.second) {
                dp[neighbour.first] = neighbour.second;
                parent[neighbour.first] = node;
                q.push(make_pair(dp[neighbour.first], neighbour.first));
            }
        }
    }
    vector<pair<int,int>> ans;
    for (int node = 1; node < numNodes; node++)
    {
        ans.push_back(make_pair(node, parent[node]));
    }
    return ans;
}
void WeightedGraph::display()
{
    for (int i = 0; i < numNodes; i++)
    {
        cout << i << " -> ";
        for (int j = 0; j < adjList[i].size(); j++)
        {
            cout << adjList[i][j].first << " " << adjList[i][j].second << "  ";
        }
        cout << endl;
    }
}

class WeightedDirectedGraph : public WeightedGraph
{
public:
    WeightedDirectedGraph(int vertices);
    void addEdge(int v, int w, int weight);
    void removeEdge(int v, int w);
};
WeightedDirectedGraph::WeightedDirectedGraph(int vertices) : WeightedGraph(vertices) {}
void WeightedDirectedGraph::addEdge(int v, int w, int weight)
{
    adjList[v].push_back(make_pair(w, weight));
}
void WeightedDirectedGraph::removeEdge(int v, int w)
{
    for (auto itr = adjList[v].begin(); itr != adjList[v].end(); itr++)
    {
        if (itr->first == w)
        {
            adjList[v].erase(itr);
            break;
        }
    }
}

int main()
{
    WeightedGraph g(8);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(2, 3, 1);
    g.addEdge(4, 5, 1);
    g.addEdge(5, 7, 1);
    g.addEdge(4, 6, 1);
    g.addEdge(6, 7, 1);
    g.addEdge(3, 4, 1);
    
    cout << "Graph :\n";
    cout << g;

    cout << '\n';
    cout << "Number of Connected Components: ";
    cout << g.numConnectedComponents();
    
    return 0;
}