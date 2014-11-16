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
