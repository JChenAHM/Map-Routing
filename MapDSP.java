public class MapDSP {
	private double[] distTo; // distTo[v] = distance of shortest s->v path
	private DirectedEdge[] edgeTo; // edgeTo[v] = last edge on shortest s->v
									// path
	private MapsPQ<Double> pq; // priority queue of vertices

	private Stack<Integer> onPQ; // vertices in pq
	private EdgeWeightedDigraph G;
	private Point2D[] points;

	/**
	 * If s = -1 is passed, the constructor simply initializes the distTo and
	 * edgeTo arrays, ready for calls to shortestPath
	 */
	public MapDSP(EdgeWeightedDigraph G, Point2D[] points, int s) {
		// initialize instance variables
		distTo = new double[G.V()];
		edgeTo = new DirectedEdge[G.V()];
		for (int v = 0; v < G.V(); v++)
			distTo[v] = Double.POSITIVE_INFINITY;
		pq = new MapsPQ<Double>(G.V());

		onPQ = new Stack<Integer>();
		this.G = new EdgeWeightedDigraph(G);
		if (points.length != G.V())
			throw new IllegalArgumentException(
					"The number of vertices on the graph must be the same as the number of points on the map.");
		this.points = new Point2D[G.V()];
		for (int i = 0; i < points.length; i++)
			this.points[i] = new Point2D(points[i].x(), points[i].y());

		if (s == -1)
			return; // don't continue if client wants to make calls to
					// shortestPath

		distTo[s] = 0.0;
		// relax vertices in order of distance from s
		pq.insert(s, distTo[s]);
		while (!pq.isEmpty()) {
			int v = pq.delMin();
			for (DirectedEdge e : G.adj(v))
				relax(e);
		}

		// check optimality conditions
		assert check(G, s);
	}

	/**
	 * Compute shortest path from source vertex to destination vertex. This only
	 * works if the second constructor has been used. This method computes one
	 * shortest path, and stops the computation as soon as the shortest path
	 * from s to d has been found. Then it only reinitializes the values in the
	 * distTo and edgeTo arrays that have changed
	 */
	public double shortestPath(int s, int d) {
		// reset priority queue without re-instantiating
		pq.reset();
		while (!onPQ.isEmpty()) {
			int v = onPQ.pop();
			// only reset distTo and edgeTo values that have changed
			distTo[v] = Double.POSITIVE_INFINITY;
			edgeTo[v] = null;
		}

		distTo[s] = points[s].distanceTo(points[d]);
		onPQ.push(s);

		// relax vertices in order of distance from s
		pq.insert(s, distTo[s]);
		while (!pq.isEmpty()) {
			int v = pq.delMin();
			if (v == d)
				break;
			for (DirectedEdge e : G.adj(v))
				relax(e, d);
		}

		if (!hasPathTo(d))
			return Double.POSITIVE_INFINITY;
		double dist = 0;
		for (DirectedEdge e : pathTo(d))
			dist += e.weight();
		return dist;
	}

	// relax edge e and update pq if changed
	private void relax(DirectedEdge e) {
		int v = e.from(), w = e.to();
		if (distTo[w] > distTo[v] + e.weight()) {
			distTo[w] = distTo[v] + e.weight();
			edgeTo[w] = e;
			if (pq.contains(w))
				pq.decreaseKey(w, distTo[w]);
			else {
				pq.insert(w, distTo[w]);
				onPQ.push(w);
			}
		}
	}

	/**
	 * Relax edge e and update pq if changed. Implements the A* algorithm for
	 * shortest paths on Euclidean maps. Vertices closer to the destination have
	 * lower weights than vertices further from the destination. This directs
	 * graph search towards destination vertex while maintaining correctness.
	 * This algorithm can get ruined by having two points in the same place on
	 * the graph, so we don't allow edges in the edgeTo array to point back and
	 * forth
	 */
	private void relax(DirectedEdge e, int d) {
		int v = e.from(), w = e.to();
		if (edgeTo[v] != null)
			if (edgeTo[v].from() == w)
				return;
		if (distTo[w] > distTo[v] + e.weight()
				+ points[w].distanceTo(points[d])
				- points[v].distanceTo(points[d])) {
			distTo[w] = distTo[v] + e.weight()
					+ points[w].distanceTo(points[d])
					- points[v].distanceTo(points[d]);
			edgeTo[w] = e;
			if (pq.contains(w))
				pq.decreaseKey(w, distTo[w]);
			else {
				pq.insert(w, distTo[w]);
				onPQ.push(w);
			}
		}
	}

	// is there a path from s to v?
	public boolean hasPathTo(int v) {
		return distTo[v] < Double.POSITIVE_INFINITY;
	}

	// shortest path from s to v as an Iterable, null if no such path
	public Iterable<DirectedEdge> pathTo(int v) {
		if (!hasPathTo(v))
			return null;
		Stack<DirectedEdge> path = new Stack<DirectedEdge>();

		int vertex = v;
		while (true) {
			DirectedEdge e = edgeTo[vertex];
			if (e == null)
				break;
			path.push(e);
			vertex = e.from();
		}
		return path;
	}

	// check optimality conditions:
	// (i) for all edges e: distTo[e.to()] <= distTo[e.from()] + e.weight()
	// (ii) for all edge e on the SPT: distTo[e.to()] == distTo[e.from()] +
	// e.weight()
	private boolean check(EdgeWeightedDigraph G, int s) {

		// check that edge weights are nonnegative
		for (DirectedEdge e : G.edges()) {
			if (e.weight() < 0) {
				System.err.println("negative edge weight detected");
				return false;
			}
		}

		// check that distTo[v] and edgeTo[v] are consistent
		if (distTo[s] != 0.0 || edgeTo[s] != null) {
			System.err.println("distTo[s] and edgeTo[s] inconsistent");
			return false;
		}
		for (int v = 0; v < G.V(); v++) {
			if (v == s)
				continue;
			if (edgeTo[v] == null && distTo[v] != Double.POSITIVE_INFINITY) {
				System.err.println("distTo[] and edgeTo[] inconsistent");
				return false;
			}
		}

		// check that all edges e = v->w satisfy distTo[w] <= distTo[v] +
		// e.weight()
		for (int v = 0; v < G.V(); v++) {
			for (DirectedEdge e : G.adj(v)) {
				int w = e.to();
				if (distTo[v] + e.weight() < distTo[w]) {
					System.err.println("edge " + e + " not relaxed");
					return false;
				}
			}
		}

		// check that all edges e = v->w on SPT satisfy distTo[w] == distTo[v] +
		// e.weight()
		for (int w = 0; w < G.V(); w++) {
			if (edgeTo[w] == null)
				continue;
			DirectedEdge e = edgeTo[w];
			int v = e.from();
			if (w != e.to())
				return false;
			if (distTo[v] + e.weight() != distTo[w]) {
				System.err.println("edge " + e + " on shortest path not tight");
				return false;
			}
		}
		return true;
	}
}
