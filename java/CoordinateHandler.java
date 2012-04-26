public class CoordinateHandler {

	private Double axis;
	private Double flattening;
	private Double centralMeridian;
	private Double latOfOrigin;
	private Double scale;
	private Double falseNorthing;
	private Double falseEasting;
	
	public CoordinateHandler() {
		grs80Params();
		rt90Projection();
	}
	
	private void grs80Params() {
		axis = 6378137.0; // GRS 80.
	    flattening = 1.0 / 298.257222101; // GRS 80.
	    centralMeridian = null;
	    latOfOrigin = 0.0;
	}
	
	public void rt90Projection() {
		centralMeridian = 15.0 + 48.0 / 60.0 + 22.624306 / 3600.0;
	    scale = 1.00000561024;
	    falseNorthing = -667.711;
	    falseEasting = 1500064.274;
	}
	
	private double mathSinh(double value) {
		return 0.5 * (Math.exp(value) - Math.exp(-value));
	}
	
	private double mathCosh(double value) {
		return 0.5 * (Math.exp(value) + Math.exp(-value));
	}
	
	private double mathAtanh(double value) {
		return 0.5 * Math.log((1.0 + value) / (1.0 - value));
	}
	
	public double[] geodeticToGrid(int lat, int lon) {
		double[] x_y = new double[2];
		if(centralMeridian == null) {
			return x_y;
		}
		
		double e2 = flattening * (2.0 - flattening),
	              n = flattening / (2.0 - flattening),
	              a_roof = axis / (1.0 + n) * (1.0 + n * n / 4.0 + n * n * n * n / 64.0),
	              A = e2,
	              B = (5.0 * e2 * e2 - e2 * e2 * e2) / 6.0,
	              C = (104.0 * e2 * e2 * e2 - 45.0 * e2 * e2 * e2 * e2) / 120.0,
	              D = (1237.0 * e2 * e2 * e2 * e2) / 1260.0,
	              beta1 = n / 2.0 - 2.0 * n * n / 3.0 + 5.0 * n * n * n / 16.0 + 41.0 * n * n * n * n / 180.0,
	              beta2 = 13.0 * n * n / 48.0 - 3.0 * n * n * n / 5.0 + 557.0 * n * n * n * n / 1440.0,
	              beta3 = 61.0 * n * n * n / 240.0 - 103.0 * n * n * n * n / 140.0,
	              beta4 = 49561.0 * n * n * n * n / 161280.0,
	              
	              deg_to_rad = Math.PI / 180.0,
	              phi = lat * deg_to_rad,
	              lambda = lon * deg_to_rad,
	              lambda_zero = centralMeridian * deg_to_rad,

	              phi_star = phi - Math.sin(phi) * Math.cos(phi) * (A +
	                      B * Math.pow(Math.sin(phi), 2) +
	                      C * Math.pow(Math.sin(phi), 4) +
	                      D * Math.pow(Math.sin(phi), 6)),
	              delta_lambda = lambda - lambda_zero,
	              xi_prim = Math.atan(Math.tan(phi_star) / Math.cos(delta_lambda)),
	              eta_prim = mathAtanh(Math.cos(phi_star) * Math.sin(delta_lambda)),
	              x = scale * a_roof * (xi_prim +
	                      beta1 * Math.sin(2.0 * xi_prim) * mathCosh(2.0 * eta_prim) +
	                      beta2 * Math.sin(4.0 * xi_prim) * mathCosh(4.0 * eta_prim) +
	                      beta3 * Math.sin(6.0 * xi_prim) * mathCosh(6.0 * eta_prim) +
	                      beta4 * Math.sin(8.0 * xi_prim) * mathCosh(8.0 * eta_prim)) +
	                      falseNorthing,
	              y = scale * a_roof * (eta_prim +
	                      beta1 * Math.cos(2.0 * xi_prim) * mathSinh(2.0 * eta_prim) +
	                      beta2 * Math.cos(4.0 * xi_prim) * mathSinh(4.0 * eta_prim) +
	                      beta3 * Math.cos(6.0 * xi_prim) * mathSinh(6.0 * eta_prim) +
	                      beta4 * Math.cos(8.0 * xi_prim) * mathSinh(8.0 * eta_prim)) +
	                      falseEasting;
		
		x_y[0] = Math.round(x * 1000.0) / 1000.0;
		x_y[1] = Math.round(y * 1000.0) / 1000.0;
		
		return x_y;
	}
	
	public double[] gridToGeodetic(int x, int y) {
		double[] lat_lon = new double[2];
		if (centralMeridian == null) {
		  return lat_lon;
		}
		
		// Prepare ellipsoid-based stuff.
		double e2 = flattening * (2.0 - flattening);
		double n = flattening / (2.0 - flattening);
		double a_roof = axis / (1.0 + n) * (1.0 + n * n / 4.0 + n * n * n * n / 64.0);
		double delta1 = n / 2.0 - 2.0 * n * n / 3.0 + 37.0 * n * n * n / 96.0 - n * n * n * n / 360.0;
		double delta2 = n * n / 48.0 + n * n * n / 15.0 - 437.0 * n * n * n * n / 1440.0;
		double delta3 = 17.0 * n * n * n / 480.0 - 37 * n * n * n * n / 840.0;
		double delta4 = 4397.0 * n * n * n * n / 161280.0;
		
		double Astar = e2 + e2 * e2 + e2 * e2 * e2 + e2 * e2 * e2 * e2;
		double Bstar = -(7.0 * e2 * e2 + 17.0 * e2 * e2 * e2 + 30.0 * e2 * e2 * e2 * e2) / 6.0;
		double Cstar = (224.0 * e2 * e2 * e2 + 889.0 * e2 * e2 * e2 * e2) / 120.0;
		double Dstar = -(4279.0 * e2 * e2 * e2 * e2) / 1260.0;
		
		// Convert
		double deg_to_rad = Math.PI / 180;
		double lambda_zero = centralMeridian * deg_to_rad;
		double xi = (x - falseNorthing) / (scale * a_roof);
		double eta = (y - falseEasting) / (scale * a_roof);
		double xi_prim = xi -
				delta1 * Math.sin(2.0 * xi) * mathCosh(2.0 * eta) -
				delta2 * Math.sin(4.0 * xi) * mathCosh(4.0 * eta) -
				delta3 * Math.sin(6.0 * xi) * mathCosh(6.0 * eta) -
				delta4 * Math.sin(8.0 * xi) * mathCosh(8.0 * eta);
		double eta_prim = eta -
				delta1 * Math.cos(2.0 * xi) * mathSinh(2.0 * eta) -
				delta2 * Math.cos(4.0 * xi) * mathSinh(4.0 * eta) -
				delta3 * Math.cos(6.0 * xi) * mathSinh(6.0 * eta) -
				delta4 * Math.cos(8.0 * xi) * mathSinh(8.0 * eta);
		double phi_star = Math.asin(Math.sin(xi_prim) / mathCosh(eta_prim));
		double delta_lambda = Math.atan(mathSinh(eta_prim) / Math.cos(xi_prim));
		double lon_radian = lambda_zero + delta_lambda;
		double lat_radian = phi_star + Math.sin(phi_star) * Math.cos(phi_star) *
				(Astar +
						Bstar * Math.pow(Math.sin(phi_star), 2) +
						Cstar * Math.pow(Math.sin(phi_star), 4) +
						Dstar * Math.pow(Math.sin(phi_star), 6));
		
		lat_lon[0] = lat_radian * 180.0 / Math.PI;
		lat_lon[1] = lon_radian * 180.0 / Math.PI;
		
		return lat_lon;
	}
}
