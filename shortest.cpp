#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <random>
using namespace std;

// 定义节点数量
#define V 36

// 定义边的结构体
struct Edge {
    int src, dest, weight;
};

int generateRandomNumber() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(1, 4);
    return distribution(gen);
}

// 定义 FatTree 类
class FatTree {
public:
    int k; // 参数 k
    int coreSwitches; // 核心层交换机数量
    int aggregationSwitches; // 汇聚层交换机数量
    int edgeSwitches; // 接入层交换机数量
    int serversPerSwitch; // 每个接入层交换机连接的服务器数量

    // 构造函数
    FatTree(int k) : k(k) {
        coreSwitches = k * k / 4;         //4
        aggregationSwitches = k * k / 2;  //8
        edgeSwitches = k * k / 2;         //8
        serversPerSwitch = 2;
    }

    // 生成邻接矩阵
    vector<vector<int>> generateAdjacencyMatrix() {
        int totalSwitches = coreSwitches + aggregationSwitches + edgeSwitches + edgeSwitches * serversPerSwitch;
        vector<vector<int>> adjacencyMatrix(totalSwitches, vector<int>(totalSwitches, 0));

        // 连接核心层和汇聚层
        int aggregationStartIndex = coreSwitches;
        for (int i = 0; i < coreSwitches; ++i) {
            for (int j = 0; j < k; ++j) {
                adjacencyMatrix[i][aggregationStartIndex + j * k / 2 + (i / (k / 2))] = generateRandomNumber();
                adjacencyMatrix[aggregationStartIndex + j * k / 2 + (i / (k / 2))][i] = generateRandomNumber();
            }
        }

        // 连接汇聚层和接入层，并连接接入层与服务器
        int edgeStartIndex = coreSwitches + aggregationSwitches;
        int serverStartIndex = coreSwitches + aggregationSwitches + edgeSwitches;
        for (int i = aggregationStartIndex; i < edgeStartIndex; ++i) {
            for (int l = 0; l < k / 2; ++l) {
                int edgeIndex = 2 * (i / 2) + aggregationSwitches + l;
                adjacencyMatrix[i][edgeIndex] = generateRandomNumber();
                adjacencyMatrix[edgeIndex][i] = generateRandomNumber();
            }
            // 连接每个接入层交换机和2台服务器
            for (int m = 0; m < serversPerSwitch; ++m) {
                int serverIndex = edgeStartIndex + 2 * i + m;
                adjacencyMatrix[i + edgeSwitches][serverIndex] = generateRandomNumber();
                adjacencyMatrix[serverIndex][i + edgeSwitches] = generateRandomNumber();
            }

        }

        return adjacencyMatrix;
    }
};

// 找到距离值数组中最小距离的节点
int minDistance(vector<int>& dist, vector<bool>& sptSet) {
    int minDist = INT_MAX, minIndex;

    for (int v = 0; v < V; v++) {
        if (!sptSet[v] && dist[v] <= minDist) {
            minDist = dist[v];
            minIndex = v;
        }
    }

    return minIndex;
}

// 打印最终的解决方案和最短路径
void printSolution(vector<int>& dist, vector<int>& parent, int src) {
    cout << "节点到源节点的最短距离和路径：" << endl;
    for (int i = 20; i < V; i++) {
        if (i != src) {
            cout << "从 " << src << " 到 " << i << " 的最短距离为：" << dist[i] << "，路径为： ";
            int current = i;
            while (current != src) {
                cout << current << " <- ";
                current = parent[current];
            }
            cout << src << endl;
        }
    }
}

// 检测负权回路
bool hasNegativeCycle(vector<Edge>& edges, vector<int>& dist) {
    for (const Edge& edge : edges) {
        int u = edge.src;
        int v = edge.dest;
        int weight = edge.weight;

        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            return true;  // 存在负权回路
        }
    }
    return false;  // 不存在负权回路
}

// 执行Dijkstra算法
void dijkstra(vector<vector<int>>& graph, int src) {
    vector<int> dist(V, INT_MAX);  // 存储最短距离的数组
    vector<int> parent(V, -1);  // 存储最短路径的父节点
    vector<bool> sptSet(V, false);  // 记录节点是否在最短路径树中

    dist[src] = 0;  // 源节点到自身的距离为0

    for (int count = 0; count < V - 1; count++) {
        int u = minDistance(dist, sptSet);  // 选择距离值最小的节点

        sptSet[u] = true;  // 将选定的节点标记为已处理

        // 更新相邻节点的最短距离和最短路径父节点
        for (int v = 0; v < V; v++) {
            if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX &&
                dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
                parent[v] = u;  // 更新最短路径的父节点
            }
        }
    }

    // 打印最终解决方案和最短路径
    printSolution(dist, parent, src);
}

// 执行Bellman-Ford算法
void bellmanFord(vector<Edge>& edges, int src) {
    vector<int> dist(V, INT_MAX);  // 存储最短距离的数组
    vector<int> parent(V, -1);  // 存储最短路径的父节点
    dist[src] = 0;  // 源节点到自身的距离为0

    // 松弛操作，对每条边进行V-1次松弛操作
    for (int i = 0; i < V - 1; i++) {
        for (const Edge& edge : edges) {
            int u = edge.src;
            int v = edge.dest;
            int weight = edge.weight;

            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                parent[v] = u;  // 更新最短路径的父节点
            }
        }
    }

    // 打印最终解决方案和最短路径
    if (hasNegativeCycle(edges, dist)) {
        cout << "图中存在负权回路，Bellman-Ford算法无法处理。" << endl;
    }
    else {
        printSolution(dist, parent, src);
    }
}

int main() {
    // 生成随机图数据（示例中使用邻接矩阵表示）
    int k = 4;
    FatTree fatTree(k);

    vector<vector<int>> graph = fatTree.generateAdjacencyMatrix();

    // 输出邻接矩阵
    for (const auto& row : graph) {
        for (int value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    // 生成边数据
    vector<Edge> edges;
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            if (graph[i][j] != 0) {
                edges.push_back({ i, j, graph[i][j] });
                edges.push_back({ j, i, graph[j][i] });
            }
        }
    }

    // 设置源节点
    int src = 20;

    // 运行Dijkstra算法，并记录执行时间
    auto start_time_dijkstra = chrono::high_resolution_clock::now();
    dijkstra(graph, src);
    auto end_time_dijkstra = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time_dijkstra = end_time_dijkstra - start_time_dijkstra;
    cout << "Dijkstra算法执行时间：" << elapsed_time_dijkstra.count() << " 秒" << endl;

    // 运行Bellman-Ford算法，并记录执行时间
    auto start_time_bellman_ford = chrono::high_resolution_clock::now();
    bellmanFord(edges, src);
    auto end_time_bellman_ford = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time_bellman_ford = end_time_bellman_ford - start_time_bellman_ford;
    cout << "Bellman-Ford算法执行时间：" << elapsed_time_bellman_ford.count() << " 秒" << endl;


    return 0;
}
