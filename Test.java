public class Test {

	public static void main(String[] args) {
		MapDSP dsp;
		EdgeWeightedDigraph g;
		int V;
		int trials;

		In in = new In(args[0]);
		trials = Integer.parseInt(args[1]);

		V = Integer.parseInt(in.readLine().trim().split("\\s+")[0]);
		Point2D[] places = new Point2D[V];
		g = new EdgeWeightedDigraph(V);

		int p1;
		int p2;
		double dist;
		while (!in.isEmpty()) {
			String[] next = in.readLine().trim().split("\\s+");
			// split with white space as delimiter

			if (next.length == 3) {
				places[Integer.parseInt(next[0])] = new Point2D(
						Double.parseDouble(next[1]),
						Double.parseDouble(next[2]));
			}

			if (next.length == 2) {
				p1 = Integer.parseInt(next[0]);
				p2 = Integer.parseInt(next[1]);
				dist = places[p1].distanceTo(places[p2]);
				g.addEdge(new DirectedEdge(p1, p2, dist));
				g.addEdge(new DirectedEdge(p2, p1, dist));
			}
		}

		dsp = new MapDSP(g, places, -1);
		for (int i = 0; i < trials; i++) {
			int s = (int) (Math.random() * V);
			int d = (int) (Math.random() * V);
			StdOut.println("from " + s + " to " + d + " : "
					+ dsp.shortestPath(s, d));

			if (dsp.hasPathTo(d)) {
				for (DirectedEdge e : dsp.pathTo(d))
					StdOut.print(e.from() + " --> ");
				StdOut.println(d);
			}
		}
	}
}
