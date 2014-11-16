public class Test {

	public static void main(String[] args) {
		EdgeWeightedDigraph g;
		int V = 0;
		int E = 0;
		int count = 0;

		Readfile in = new Readfile("input.txt");
		
        String[] next1 = in.readLine().trim().split("\\s+");
		V = Integer.parseInt(next1[0]);
		E = Integer.parseInt(next1[1]);
		
		Point2D[] places = new Point2D[V];
		g = new EdgeWeightedDigraph(V);

		int p1;
		int p2;
		int s;
		int d;
		double dist;
		
		while (!in.isEmpty()) {	
			String[] next2 = in.readLine().trim().split("\\s+");
			if (next2.length == 3) {
				places[Integer.parseInt(next2[0])] = new Point2D(
						Double.parseDouble(next2[1]),
						Double.parseDouble(next2[2]));
			}

			if (next2.length == 2 && count < E) {
				count++;
				p1 = Integer.parseInt(next2[0]);
				p2 = Integer.parseInt(next2[1]);
				dist = places[p1].distanceTo(places[p2]);
				g.addEdge(new DirectedEdge(p1, p2, dist));
				g.addEdge(new DirectedEdge(p2, p1, dist));
			}
			
			if (next2.length == 2 && count >= E){
				Dijkstra dsp = new Dijkstra(g, places, -1);
				s = Integer.parseInt(next2[0]);
				d = Integer.parseInt(next2[1]);
				System.out.println(dsp.shortestPath(s, d));
				if (dsp.hasPathTo(d)) {
					System.out.print("[");
					for (DirectedEdge e : dsp.pathTo(d))
						System.out.print(e.from() + ",");
					System.out.println(d+ "]");
				}
				count++;
			}
		}

		/*dsp = new Dijkstra(g, places, -1);
		int s1 = 0;
		int d1 = 5;
		int s2 = 2;
		int d2 = 3;
		System.out.println(dsp.shortestPath(s1, d1));
		if (dsp.hasPathTo(d1)) {
			System.out.print("[");
			for (DirectedEdge e : dsp.pathTo(d1))
				System.out.print(e.from() + ",");
			System.out.println(d1+ "]");
		}
		System.out.println(dsp.shortestPath(s2, d2));
		if (dsp.hasPathTo(d2)) {
			System.out.print("[");
			for (DirectedEdge e : dsp.pathTo(d2))
				System.out.print(e.from() + ",");
			System.out.println(d2+ "]");
		}*/
	}
}
