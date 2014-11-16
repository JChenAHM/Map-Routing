import java.util.Iterator;
import java.util.NoSuchElementException;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.net.URL;
import java.net.URLConnection;
import java.util.Locale;
import java.util.Scanner;

public class Dijkstra {
	
	/**
	 * Bag
	 */
	public static class Bag<Item> implements Iterable<Item> {
	    private int N;               // number of elements in bag
	    private Node<Item> first;    // beginning of bag

	    /*
	     * helper linked list class
	     */
	    @SuppressWarnings("hiding")
		private class Node<Item> {
	        private Item item;
	        private Node<Item> next;
	    }
	  
	    /*
	     * Adds the item to this bag.
	     */
	    public void add(Item item) {
	        Node<Item> oldfirst = first;
	        first = new Node<Item>();
	        first.item = item;
	        first.next = oldfirst;
	        N++;
	    }

	    /*
	     * Returns an iterator that iterates over the items in the bag in arbitrary order.
	     */
	    public Iterator<Item> iterator()  {
	        return new ListIterator<Item>(first);  
	    }

	    @SuppressWarnings("hiding")
		private class ListIterator<Item> implements Iterator<Item> {
	        private Node<Item> current;

	        public ListIterator(Node<Item> first) {
	            current = first;
	        }

	        public boolean hasNext()  { return current != null;                     }
	        public void remove()      { throw new UnsupportedOperationException();  }

	        public Item next() {
	            if (!hasNext()) throw new NoSuchElementException();
	            Item item = current.item;
	            current = current.next; 
	            return item;
	        }
	    }
	}
	
	/**
	 * DirectedEdge
	 */
	public static class DirectedEdge { 
	    private final int v;
	    private final int w;
	    private final double weight;

	    /*
	     * Initializes a directed edge from vertex v to vertex w with the given weight.
	     */
	    public DirectedEdge(int v, int w, double weight) {
	        if (v < 0) throw new IndexOutOfBoundsException("Vertex names must be nonnegative integers");
	        if (w < 0) throw new IndexOutOfBoundsException("Vertex names must be nonnegative integers");
	        if (Double.isNaN(weight)) throw new IllegalArgumentException("Weight is NaN");
	        this.v = v;
	        this.w = w;
	        this.weight = weight;
	    }

	    /*
	     * Returns the tail vertex of the directed edge.
	     */
	    public int from() {
	        return v;
	    }

	    /*
	     * Returns the head vertex of the directed edge.
	     */
	    public int to() {
	        return w;
	    }

	    /*
	     * Returns the weight of the directed edge.
	     */
	    public double weight() {
	        return weight;
	    }
	}
	
	/**
	 * Stack
	 */
	public static class Stack<Item> implements Iterable<Item> {
	    private int N;                // size of the stack
	    private Node<Item> first;     // top of stack

	    // helper linked list class
	    @SuppressWarnings("hiding")
		private class Node<Item> {
	        private Item item;
	        private Node<Item> next;
	    }
	    
	    public boolean isEmpty() {
	        return first == null;
	    }   

	    /*
	     * Adds the item to this stack.
	     */
	    public void push(Item item) {
	        Node<Item> oldfirst = first;
	        first = new Node<Item>();
	        first.item = item;
	        first.next = oldfirst;
	        N++;
	    }
	    
	    /*
	     * Removes and returns the item most recently added to this stack
	     */
	    public Item pop() {
	        if (isEmpty()) throw new NoSuchElementException("Stack underflow");
	        Item item = first.item;        // save item to return
	        first = first.next;            // delete first node
	        N--;
	        return item;                   // return the saved item
	    }
	 
	    /*
	     * Returns an iterator to this stack that iterates through the items in LIFO order.
	     */
	    public Iterator<Item> iterator() {
	        return new ListIterator<Item>(first);
	    }

	    @SuppressWarnings("hiding")
		private class ListIterator<Item> implements Iterator<Item> {
	        private Node<Item> current;

	        public ListIterator(Node<Item> first) {
	            current = first;
	        }
	        public boolean hasNext()  { return current != null;                     }
	        public void remove()      { throw new UnsupportedOperationException();  }

	        public Item next() {
	            if (!hasNext()) throw new NoSuchElementException();
	            Item item = current.item;
	            current = current.next; 
	            return item;
	        }
	    }
	}

	/**
	 * EdgeWeightedDigraph
	 */
	public static class EdgeWeightedDigraph {
	    private final int V;
	    private int E;
	    private Bag<DirectedEdge>[] adj;
	    
	    /*
	     * Create an empty edge-weighted digraph with V vertices.
	     */
	    @SuppressWarnings("unchecked")
		public EdgeWeightedDigraph(int V) {
	        if (V < 0) throw new IllegalArgumentException("Number of vertices in a Digraph must be nonnegative");
	        this.V = V;
	        this.E = 0;
	        adj = (Bag<DirectedEdge>[]) new Bag[V];
	        for (int v = 0; v < V; v++)
	            adj[v] = new Bag<DirectedEdge>();
	    }

	    /*
	     * Copy constructor.
	     */
	    public EdgeWeightedDigraph(EdgeWeightedDigraph G) {
	        this(G.V());
	        this.E = G.E();
	        for (int v = 0; v < G.V(); v++) {
	            // reverse so that adjacency list is in same order as original
	            Stack<DirectedEdge> reverse = new Stack<DirectedEdge>();
	            for (DirectedEdge e : G.adj[v]) {
	                reverse.push(e);
	            }
	            for (DirectedEdge e : reverse) {
	                adj[v].add(e);
	            }
	        }
	    }

	    /*
	     * Return the number of vertices in this digraph.
	     */
	    public int V() {
	        return V;
	    }

	    /*
	     * Return the number of edges in this digraph.
	     */
	    public int E() {
	        return E;
	    }


	    /*
	     * Add the directed edge e to this digraph.
	     */
	    public void addEdge(DirectedEdge e) {
	        int v = e.from();
	        adj[v].add(e);
	        E++;
	    }

	    public Iterable<DirectedEdge> adj(int v) {
	        return adj[v];
	    }

	    /*
	     * Return all edges in this digraph as an Iterable.
	     */
	    public Iterable<DirectedEdge> edges() {
	        Bag<DirectedEdge> list = new Bag<DirectedEdge>();
	        for (int v = 0; v < V; v++) {
	            for (DirectedEdge e : adj(v)) {
	                list.add(e);
	            }
	        }
	        return list;
	    } 
	}
	
	/**
	 * MapsPQ
	 */
	public class MapsPQ<Key extends Comparable<Key>> implements Iterable<Integer> {
        private int NMAX; // maximum number of elements on PQ
        private int N; // number of elements on PQ
        private int[] pq; // binary heap using 1-based indexing
        private int[] qp; // inverse of pq -- qp[pq[i]] = pq[qp[i]] = i
        private Key[] keys; // keys[i] = priority of i

        /*
         * Create an empty indexed priority queue with indices between 0 and NMAX-1.
         */
        @SuppressWarnings("unchecked")
		public MapsPQ(int NMAX) {
        	if (NMAX < 0)
        		throw new IllegalArgumentException();
        	this.NMAX = NMAX;
        	keys = (Key[]) new Comparable[NMAX + 1]; // make this of length NMAX??
        	pq = new int[NMAX + 1];
        	qp = new int[NMAX + 1]; // make this of length NMAX??
        	for (int i = 0; i <= NMAX; i++)
        		qp[i] = -1;
        }

        public boolean isEmpty() {
        	return N == 0;
        } 

        /*
         * Is i an index on the priority queue?
         */
        public boolean contains(int i) {
        	if (i < 0 || i >= NMAX)
        		throw new IndexOutOfBoundsException();
        	return qp[i] != -1;
        }

        /*
         * Associate key with index i.
         */
        public void insert(int i, Key key) {
        	if (i < 0 || i >= NMAX)
        		throw new IndexOutOfBoundsException();
        	if (contains(i))
        		throw new IllegalArgumentException("index is already in the priority queue");
        	N++;
        	qp[i] = N;
        	pq[N] = i;
        	keys[i] = key;
        	swim(N);
        }

        /*
         * Reset PQ. If the queue is not full, this runs faster than instantiating a new PQ
         */
        public void reset() {
        	int v;
        	for (int i = N; i > 0; i--) {
        		v = pq[i];
        		pq[i] = 0;
        		qp[v] = -1;
        		keys[v] = null;
        		// resetting keys not necessary
        	}
        	
        	N = 0;
        }

        /*
         * Delete a minimal key and return its associated index.
         */
        public int delMin() {
        	if (N == 0)
        		throw new NoSuchElementException("Priority queue underflow");
        	int min = pq[1];
        	exch(1, N--);
        	sink(1);
        	qp[min] = -1; // delete
        	keys[pq[N + 1]] = null; // to help with garbage collection
        	pq[N + 1] = 0; // delete
        	return min;
        }

        /*
         * Decrease the key associated with index i to the specified value.
         */
        public void decreaseKey(int i, Key key) {
        	if (i < 0 || i >= NMAX)
        		throw new IndexOutOfBoundsException();
        	if (!contains(i))
        		throw new NoSuchElementException("index is not in the priority queue");
        	if (keys[i].compareTo(key) <= 0)
        		throw new IllegalArgumentException("Calling decreaseKey() with given argument would not strictly decrease the key");
        	keys[i] = key;
        	swim(qp[i]);
        }

        private boolean greater(int i, int j) {
        	return keys[pq[i]].compareTo(keys[pq[j]]) > 0;
        }
        
        private void exch(int i, int j) {
        	int swap = pq[i];
        	pq[i] = pq[j];
        	pq[j] = swap;
        	qp[pq[i]] = i;
        	qp[pq[j]] = j;
        }

        private void swim(int k) {
        	while (k > 1 && greater(k / 2, k)) {
        		exch(k, k / 2);
        		k = k / 2;
        	}
        }

        private void sink(int k) {
        	while (2 * k <= N) {
        		int j = 2 * k;
        		if (j < N && greater(j, j + 1))
        			j++;
        		if (!greater(k, j))
        			break;
        		exch(k, j);
        		k = j;
        	}
        }

        /*
         * Return an iterator that iterates over all of the elements on the priority queue in ascending order.
         */
        public Iterator<Integer> iterator() {
	        return new HeapIterator();
        }

        private class HeapIterator implements Iterator<Integer> {
        	// create a new pq
        	private MapsPQ<Key> copy;

        	// add all elements to copy of heap
        	// takes linear time since already in heap order so no keys move
        	public HeapIterator() {
        		copy = new MapsPQ<Key>(pq.length - 1);
        		for (int i = 1; i <= N; i++)
        			copy.insert(pq[i], keys[pq[i]]);
        	}
        	
        	public boolean hasNext() {
        		return !copy.isEmpty();
        	}
        	
        	public void remove() {
        		throw new UnsupportedOperationException();
        	}

            public Integer next() {
            	if (!hasNext())
            		throw new NoSuchElementException();
            	return copy.delMin();
            }
        }
    }
	
	/**
	 * Point2D
	 */
	public static class Point2D {

	    private final double x;    // x coordinate
	    private final double y;    // y coordinate

	    // create a new point (x, y)
	    public Point2D(double x, double y) {
	        this.x = x;
	        this.y = y;
	    }

	    // return the x-coorindate of this point
	    public double x() { return x; }

	    // return the y-coorindate of this point
	    public double y() { return y; }

	    // return Euclidean distance between this point and that point
	    public double distanceTo(Point2D that) {
	        double dx = this.x - that.x;
	        double dy = this.y - that.y;
	        double distance = Math.sqrt(dx*dx + dy*dy);
	        BigDecimal bg = new BigDecimal(distance);
	        double distance1 = bg.setScale(2, BigDecimal.ROUND_HALF_UP).doubleValue();
	        return distance1;
	    }
	}

	/**
	 * Dijkstra
	 */
	private double[] distTo; // distTo[v] = distance of shortest s->v path
	private DirectedEdge[] edgeTo; // edgeTo[v] = last edge on shortest s->v path
	private MapsPQ<Double> pq; // priority queue of vertices

	private Stack<Integer> onPQ; // vertices in pq
	private EdgeWeightedDigraph G;
	private Point2D[] points;

	/*
	 * If s = -1 is passed, the constructor simply initializes the distTo and
	 * edgeTo arrays, ready for calls to shortestPath
	 */
	public Dijkstra(EdgeWeightedDigraph G, Point2D[] points, int s) {
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
			return; // don't continue if client wants to make calls to shortestPath

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

	/*
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

	/*
	 *  relax edge e and update pq if changed
	 */
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

	/*
	 * Relax edge e and update pq if changed. Implements the A* algorithm for
	 * shortest paths on Euclidean maps. Vertices closer to the destination have
	 * lower weights than vertices further from the destination. This directs
	 * graph search towards destination vertex while maintaining correctness.
	 * This algorithm can get ruined by having two points in the same place on
	 * the graph, so we don't allow edges in the edgeTo array to point back and forth
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

	/*
	 * is there a path from s to v?
	 */
	public boolean hasPathTo(int v) {
		return distTo[v] < Double.POSITIVE_INFINITY;
	}

	/*
	 * shortest path from s to v as an Iterable, null if no such path
	 */
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

	/*
	 * check optimality conditions:
	 * (i) for all edges e: distTo[e.to()] <= distTo[e.from()] + e.weight()
	 * (ii) for all edge e on the SPT: distTo[e.to()] == distTo[e.from()] + e.weight()
	 */
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

		// check that all edges e = v->w satisfy distTo[w] <= distTo[v] + e.weight()
		for (int v = 0; v < G.V(); v++) {
			for (DirectedEdge e : G.adj(v)) {
				int w = e.to();
				if (distTo[v] + e.weight() < distTo[w]) {
					System.err.println("edge " + e + " not relaxed");
					return false;
				}
			}
		}

		// check that all edges e = v->w on SPT satisfy distTo[w] == distTo[v] + e.weight()
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
	
    /**
     * Readfile
     */
	public static final class Readfile {
	    
	    private Scanner scanner;
	    
	    // assume Unicode UTF-8 encoding
	    private static final String CHARSET_NAME = "UTF-8";

	    // assume language = English, country = US for consistency with System.out.
	    private final Locale LOCALE = Locale.US;
	    
	    /*
	     * Create an input stream from a socket.
	     */
	    public Readfile(java.net.Socket socket) {
	        try {
	            InputStream is = socket.getInputStream();
	            scanner = new Scanner(new BufferedInputStream(is), CHARSET_NAME);
	            scanner.useLocale(LOCALE);
	        }
	        catch (IOException ioe) {
	            System.err.println("Could not open " + socket);
	        }
	    }

	    /*
	     * Create an input stream from a URL.
	     */
	    public Readfile(URL url) {
	        try {
	            URLConnection site = url.openConnection();
	            InputStream is     = site.getInputStream();
	            scanner            = new Scanner(new BufferedInputStream(is), CHARSET_NAME);
	            scanner.useLocale(LOCALE);
	        }
	        catch (IOException ioe) {
	            System.err.println("Could not open " + url);
	        }
	    }
	    
	    /*
	     * Create an input stream from a filename or web page name.
	     */
	    public Readfile(String s) {
	        try {
	            // first try to read file from local file system
	            File file = new File(s);
	            if (file.exists()) {
	                scanner = new Scanner(file, CHARSET_NAME);
	                scanner.useLocale(LOCALE);
	                return;
	            }

	            // next try for files included in jar
	            URL url = getClass().getResource(s);

	            // or URL from web
	            if (url == null) { url = new URL(s); }

	            URLConnection site = url.openConnection();

	            InputStream is     = site.getInputStream();
	            scanner            = new Scanner(new BufferedInputStream(is), CHARSET_NAME);
	            scanner.useLocale(LOCALE);
	        }
	        catch (IOException ioe) {
	            System.err.println("Could not open " + s);
	        }
	    } 
	    
	    public boolean isEmpty() {
	        return !scanner.hasNext();
	    }

	    /*
	     * Read and return the next line.
	     */
	    public String readLine() {
	        String line;
	        try                 { line = scanner.nextLine(); }
	        catch (Exception e) { line = null;               }
	        return line;
	    }
	    
	    /*
	     * Close the input stream.
	     */
	    public void close() {
	        scanner.close();  
	    }
	}
	
	/**
	 * Test
	 */
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
