var CoordinateHandler = function() {
  var axis = null,
          flattening = null,
          centralMeridian = null,
          latOfOrigin = null,
          scale = null,
          falseNorthing = null,
          falseEasting = null;

  (function grs80Params() {
    axis = 6378137.0; // GRS 80.
    flattening = 1.0 / 298.257222101; // GRS 80.
    centralMeridian = null;
    latOfOrigin = 0.0;
  })();

  //rt90_2.5_gon_v
  (function rt90Projection() {
    centralMeridian = 15.0 + 48.0 / 60.0 + 22.624306 / 3600.0;
    scale = 1.00000561024;
    falseNorthing = -667.711;
    falseEasting = 1500064.274;
  })();

  function mathSinh(value) {
    return 0.5 * (Math.exp(value) - Math.exp(-value));
  }

  function mathCosh(value) {
    return 0.5 * (Math.exp(value) + Math.exp(-value));
  }

  function mathAtanh(value) {
    return 0.5 * Math.log((1.0 + value) / (1.0 - value));
  }

  return {
    geodeticToGrid: function(lat, long) {
      var x_y = new Array(2);
      if (centralMeridian == null) {
        return x_y;
      }
      var e2 = flattening * (2.0 - flattening),
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
              lambda = long * deg_to_rad,
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

      return {x:x_y[0],y:x_y[1]};
    },
    gridToGeodetic: function (x, y) {
      var lat_lon = new Array(2);
      if (centralMeridian == null) {
        return lat_lon;
      }
      // Prepare ellipsoid-based stuff.
      var e2 = flattening * (2.0 - flattening);
      var n = flattening / (2.0 - flattening);
      var a_roof = axis / (1.0 + n) * (1.0 + n * n / 4.0 + n * n * n * n / 64.0);
      var delta1 = n / 2.0 - 2.0 * n * n / 3.0 + 37.0 * n * n * n / 96.0 - n * n * n * n / 360.0;
      var delta2 = n * n / 48.0 + n * n * n / 15.0 - 437.0 * n * n * n * n / 1440.0;
      var delta3 = 17.0 * n * n * n / 480.0 - 37 * n * n * n * n / 840.0;
      var delta4 = 4397.0 * n * n * n * n / 161280.0;

      var Astar = e2 + e2 * e2 + e2 * e2 * e2 + e2 * e2 * e2 * e2;
      var Bstar = -(7.0 * e2 * e2 + 17.0 * e2 * e2 * e2 + 30.0 * e2 * e2 * e2 * e2) / 6.0;
      var Cstar = (224.0 * e2 * e2 * e2 + 889.0 * e2 * e2 * e2 * e2) / 120.0;
      var Dstar = -(4279.0 * e2 * e2 * e2 * e2) / 1260.0;

      // Convert.
      var deg_to_rad = Math.PI / 180;
      var lambda_zero = centralMeridian * deg_to_rad;
      var xi = (x - falseNorthing) / (scale * a_roof);
      var eta = (y - falseEasting) / (scale * a_roof);
      var xi_prim = xi -
              delta1 * Math.sin(2.0 * xi) * mathCosh(2.0 * eta) -
              delta2 * Math.sin(4.0 * xi) * mathCosh(4.0 * eta) -
              delta3 * Math.sin(6.0 * xi) * mathCosh(6.0 * eta) -
              delta4 * Math.sin(8.0 * xi) * mathCosh(8.0 * eta);
      var eta_prim = eta -
              delta1 * Math.cos(2.0 * xi) * mathSinh(2.0 * eta) -
              delta2 * Math.cos(4.0 * xi) * mathSinh(4.0 * eta) -
              delta3 * Math.cos(6.0 * xi) * mathSinh(6.0 * eta) -
              delta4 * Math.cos(8.0 * xi) * mathSinh(8.0 * eta);
      var phi_star = Math.asin(Math.sin(xi_prim) / mathCosh(eta_prim));
      var delta_lambda = Math.atan(mathSinh(eta_prim) / Math.cos(xi_prim));
      var lon_radian = lambda_zero + delta_lambda;
      var lat_radian = phi_star + Math.sin(phi_star) * Math.cos(phi_star) *
              (Astar +
                      Bstar * Math.pow(Math.sin(phi_star), 2) +
                      Cstar * Math.pow(Math.sin(phi_star), 4) +
                      Dstar * Math.pow(Math.sin(phi_star), 6));
      lat_lon[0] = lat_radian * 180.0 / Math.PI;
      lat_lon[1] = lon_radian * 180.0 / Math.PI;
      //return lat_lon;
      return {lat:lat_lon[0],lon:lat_lon[1]};
    }
  }

};