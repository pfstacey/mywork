import java.util.*;

/**
 * Library for graph analysis
 * 
 * @author Chris Bailey-Kellogg, Dartmouth CS 10, Fall 2016
 * 
 */
public class GraphLibrary {
	/**
	 * Takes a random walk from a vertex, up to a given number of steps
	 * So a 0-step path only includes start, while a 1-step path includes start and one of its out-neighbors,
	 * and a 2-step path includes start, an out-neighbor, and one of the out-neighbor's out-neighbors
	 * Stops earlier if no step can be taken (i.e., reach a vertex with no out-edge)
	 *
	 * @param g     graph to walk on
	 * @param start initial vertex (assumed to be in graph)
	 * @param steps max number of steps
	 * @return a list of vertices starting with start, each with an edge to the sequentially next in the list;
	 * null if start isn't in graph
	 */
	public static <V, E> List<V> randomWalk(Graph<V, E> g, V start, int steps) {
		List<V> path = new ArrayList<V>();
		path.add(start);
		List<V> allVertices = new ArrayList<V>();
		for(V v : g.vertices()) {
			allVertices.add(v);
		}
		walkPath(allVertices, g, start, steps, path);
		return path;
	}
	/*
	* walks a random path throughout all the vertices
	 */
	public static <V, E> void walkPath(List<V> allVertices, Graph<V, E> g, V start, int steps, List<V> path){
		if(steps > 0 && start != null){
			int numNeighbors = 0;
			for(V v: g.outNeighbors(start)){
				numNeighbors++;
			}
			int randomNumber = (int)(Math.random()*numNeighbors)+1;

			path.add(allVertices.get(randomNumber));
			walkPath(allVertices, g, allVertices.get(randomNumber), (--steps), path);
		}
	}

	/**
	 * Orders vertices in decreasing order by their in-degree
	 *
	 * @param g graph
	 * @return list of vertices sorted by in-degree, decreasing (i.e., largest at index 0)
	 */
	public static <V, E> List<V> verticesByInDegree(Graph<V, E> g) {
		// TODO: your code here
		List<V> outList = new ArrayList<V>();
		for(V v : g.vertices()){
				outList.add(v);
			}
			Comparator<V> compare = new Comparator<V>() {
				@Override
				public int compare(V o1, V o2) {
					return g.inDegree(o2) - g.inDegree(o1);
				}
			};
			outList.sort(compare);
		return outList;
	}


	/** Breadth First Search
	 *
	 * @param G -- graph to search
	 * @param start -- starting vertex
	 */
	public static <V,E> Graph<V,E> BFS(Graph<V,E> G, V start) {
		try {
			Graph<V, E> output = new AdjacencyMapGraph<V, E>();
			Map<V, V> backTrack = new HashMap<V, V>(); //initialize backTrack
			backTrack.put(start, null); //load start vertex with null parent
			Set<V> visited = new HashSet<V>(); //Set to track which vertices have already been visited
			Queue<V> queue = new LinkedList<V>(); //queue to implement BFS

			output.insertVertex(start);
			queue.add(start); //enqueue start vertex
			visited.add(start); //add start to visited Set
			while (!queue.isEmpty()) { //loop until no more vertices
				V u = queue.remove(); //dequeue
				for (V v : G.outNeighbors(u)) { //loop over out neighbors
					if (!visited.contains(v)) { //if neighbor not visited, then neighbor is discovered from this vertex
						visited.add(v); //add neighbor to visited Set
						queue.add(v); //enqueue neighbor
						output.insertVertex(v);
						output.insertDirected(v, u, null);
						backTrack.put(v, u); //save that this vertex was discovered from prior vertex

					}
				}
			}
			return output;
		}catch(Exception e){
			new Exception("Invalid starting actor");
			return null;
		}
	}

	/** getPath
	 *
	 * finds the path between a given vertex and the "root" vertex or center vertex
	 *
	 * @param tree -- graph to search
	 * @param v -- path between v and starting node of the tree
	 */
	public static <V,E> List<V> getPath(Graph<V,E> tree, V v){
		try {
			ArrayList<V> path = new ArrayList<V>();
			V current = v;
			path.add(current);
			while (tree.outDegree(current) != 0) {
				for (V i : tree.outNeighbors(current)) {
					current = i;
					path.add(current);
				}
			}
			return path;
		}catch(Exception e){
			new Exception("Invalid starting actor");
			return null;
		}
	}

	/** missingVertices
	 *
	 * finds the vertices that are in graph but not in subgraph
	 *
	 * @param graph -- graph of all vertices
	 * @param subgraph -- graph of some or all of the vertices of graph
	 */
	public static <V,E> Set<V> missingVertices(Graph<V,E> graph, Graph<V,E> subgraph){
		Set<V> set = new HashSet<V>();
		for(V v: graph.vertices()){
			if(!subgraph.hasVertex(v)){
				set.add(v);
			}
		}
		return set;
	}

	/** averageSeparation
	 *
	 * finds the average separation between the root and all other vertices
	 *
	 * @param tree -- graph to search
	 * @param root -- sets root from which to search
	 */
	public static <V,E> double averageSeparation(Graph<V,E> tree, V root){
		return recursionHelper(tree,root,0.0)/(tree.numVertices()-1);
	}

	/** recursionHelper
	 *
	 * exists to help averageSeparation recurse
	 *
	 * @param tree -- graph to search
	 * @param root -- sets root from which to search
	 * @param sum -- sum of the separation of all the vertices visited
	 */
	public static <V,E> double recursionHelper(Graph<V,E> tree, V root, double sum){
		if(tree.inDegree(root)==0) return sum;
		double total = sum;
		for(V v: tree.inNeighbors(root)){
			total += recursionHelper(tree, v, sum+1);
		}
		return total;
	}

}
