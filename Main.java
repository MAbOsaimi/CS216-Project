import java.util.*;

public class Main {
    static final int INF = Integer.MAX_VALUE;

    // Function to generate a random connected graph
    static int[][] generateRandomGraph(int V) {
        int[][] graph = new int[V][V];
        Random rand = new Random();

        for (int i = 0; i < V; i++) {
            graph[i][i] = INF;

            // Make sure each vertex has at least one edge (In order to get a connected
            // graph)
            int j = (i == V - 1) ? 0 : i + 1;
            graph[i][j] = graph[j][i] = rand.nextInt(10) + 1;

            // Randomly decide if there are edges between other vertices
            for (j = i + 2; j != i && j < V; j++) {
                boolean isEdge = rand.nextBoolean();
                if (isEdge) {
                    graph[i][j] = graph[j][i] = rand.nextInt(100) + 1; // Assign random cost (Between 1-100)
                } else {
                    graph[i][j] = graph[j][i] = INF;
                }
            }
        }

        return graph;
    }

    // Function to check whether a given set of edges forms a connected graph (Using
    // Depth-first-search)
    static boolean isConnected(int[][] graph, boolean[][] subset) {
        int V = graph.length;
        boolean[] visited = new boolean[V];
        Stack<Integer> stack = new Stack<>();
        stack.push(0);
        visited[0] = true;
        while (!stack.isEmpty()) {
            int u = stack.pop();
            for (int v = 0; v < V; v++) {
                if (subset[u][v] && !visited[v]) {
                    stack.push(v);
                    visited[v] = true;
                }
            }
        }
        for (boolean v : visited) {
            if (!v)
                return false;
        }
        return true;
    }

    // Function to calculate the total cost of a graph
    static int totalCost(int[][] graph) {
        int total = 0;
        for (int i = 0; i < graph.length; i++) {
            for (int j = i + 1; j < graph.length; j++) {
                if (graph[i][j] != INF) { // If there is an edge, add the cost to the total
                    total += graph[i][j];
                }
            }
        }
        return total;

    }

    // Function to calculate the total cost of a subset of edges
    static int totalCost(int[][] graph, boolean[][] subset) {
        int total = 0;
        for (int i = 0; i < graph.length; i++) {
            for (int j = i + 1; j < graph.length; j++) {
                if (subset[i][j])
                    total += graph[i][j];
            }
        }
        return total;
    }

    // Function to generate all possible subsets of edges, find connected ones
    // and find the one with the least cost
    static int bruteForce(int[][] graph) {
        int V = graph.length;

        // Count the number of edges in the graph
        int E = 0;
        for (int i = 0; i < V; i++) {
            for (int j = i + 1; j < V; j++) {
                if (graph[i][j] != INF)
                    E++;
            }
        }

        boolean[][] subset = new boolean[V][V];
        int minCost = INF;

        for (int i = 0; i < Math.pow(2, E); i++) {
            // Generate subset of edges (From 0 to 2^E - 1)
            int e = 0;
            for (int u = 0; u < V; u++) {
                for (int v = u + 1; v < V; v++) {
                    if (graph[u][v] != INF) {
                        subset[u][v] = subset[v][u] = Math.floor((i / Math.pow(2, e))) % 2 != 0;
                        e++;
                    } else {
                        subset[u][v] = subset[v][u] = false;
                    }
                }
            }

            // Check if the subset is connected and update the minCost
            if (isConnected(graph, subset)) {
                minCost = Math.min(minCost, totalCost(graph, subset));
            }
        }
        return minCost;
    }

    // Function to find the root of the tree for a vertex
    static int find(int[] parent, int i) {
        if (parent[i] == i)
            return i;
        return find(parent, parent[i]);
    }

    // Function to perform union of two subsets
    static void union(int[] parent, int x, int y) {
        int xroot = find(parent, x);
        int yroot = find(parent, y);
        parent[yroot] = xroot;
    }

    // Function to implement MergeSort algorithm for sorting edges
    static void mergeSort(int[][] edges, int left, int right) {
        if (left < right) {
            int mid = left + (right - left) / 2;
            mergeSort(edges, left, mid);
            mergeSort(edges, mid + 1, right);
            merge(edges, left, mid, right);
        }
    }

    // Function to merge two subarrays in the array
    static void merge(int[][] edges, int left, int mid, int right) {
        int n1 = mid - left + 1;
        int n2 = right - mid;

        int[][] L = new int[n1][3];
        int[][] R = new int[n2][3];

        System.arraycopy(edges, left, L, 0, n1);
        System.arraycopy(edges, mid + 1, R, 0, n2);

        int i = 0, j = 0;
        int k = left;
        while (i < n1 && j < n2) {
            if (L[i][2] <= R[j][2]) {
                edges[k] = L[i];
                i++;
            } else {
                edges[k] = R[j];
                j++;
            }
            k++;
        }

        while (i < n1) {
            edges[k] = L[i];
            i++;
            k++;
        }

        while (j < n2) {
            edges[k] = R[j];
            j++;
            k++;
        }
    }

    // Function to implement Kruskal's algorithm
    static int kruskal(int[][] graph) {
        int V = graph.length;

        int e = 0;
        int cost = 0;

        ArrayList<int[]> edgeList = new ArrayList<>();
        for (int i = 0; i < V; i++) {
            for (int j = i + 1; j < V; j++) {
                if (graph[i][j] != INF)
                    edgeList.add(new int[] { i, j, graph[i][j] });
            }
        }

        int[][] edges = edgeList.toArray(new int[0][]);

        mergeSort(edges, 0, edges.length - 1); // Sort all the edges in non-decreasing order of their cost using merge
                                               // sort

        int[] parent = new int[V];
        for (int v = 0; v < V; ++v)
            parent[v] = v;

        for (int[] edge : edges) {
            if (e == V - 1)
                break;

            int x = find(parent, edge[0]);
            int y = find(parent, edge[1]);

            if (x != y) {
                e++;
                cost += edge[2];
                union(parent, x, y);
            }
        }

        return cost;
    }

    // Function to test each algorithm
    public static void runTest(int maxVertices) {
        for (int V = 2; V <= maxVertices; V++) {
            int[][] graph = generateRandomGraph(V);

            System.out.println("\nNumber of vertices: " + V);
            System.out.println("Total cost: " + totalCost(graph));
            // System.out.println("Graph: " + Arrays.deepToString(graph));

            long startTime, endTime;

            int minCost;

            // Run bruteForce and measure time
            startTime = System.nanoTime();
            minCost = bruteForce(graph);
            endTime = System.nanoTime();
            long durBFNano = endTime - startTime;
            double durBFMillis = durBFNano / 1000000; // Convert nanosecond to millisecond (Nanosecond / 10^6)
            System.out.println("Minimum cost of graph by bruteForce: " + minCost);
            System.out.println("Execution time of bruteForce in nanoseconds: " + durBFNano);
            System.out.println("Execution time of bruteForce in milliseconds: " + durBFMillis);

            // Run kruskal and measure time
            startTime = System.nanoTime();
            minCost = kruskal(graph);
            endTime = System.nanoTime();
            long durKNano = endTime - startTime;
            double durKMillis = durKNano / 1000000; // Convert nanosecond to millisecond (Nanosecond / 10^6)

            System.out.println("Minimum cost of graph by kruskal: " + minCost);
            System.out.println("Execution time of kruskal in nanoseconds: " + durKNano);
            System.out.println("Execution time of kruskal in milliseconds: " + durKMillis);
            double improvement = 100.0 * Math.max(0, (durBFNano / durKNano - 1));
            System.out.printf("Percentage of improvement in time by kruskal over brute force (Nanoseconds): %.2f%%\n",
                    improvement);

        }
    }

    public static void main(String[] args) {

        int[][] graph = {
                { 0, 10, 20, 30 },
                { 10, 0, 40, 50 },
                { 20, 40, 0, 20 },
                { 30, 50, 20, 0 },
        };

        System.out.println(bruteForce(graph));
        System.out.println(kruskal(graph));

        runTest(10);
    }
}
