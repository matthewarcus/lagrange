"use strict";

var Lagrange = {};

(function() {
    // Vector functions
    var Vector = {
        dot: function(v0,v1) {
            return v0[0]*v1[0] + v0[1]*v1[1];
        },
        length: function(v) {
            return Math.sqrt(Vector.dot(v, v));
        },
        add: function(v0,v1) {
            return [v0[0]+v1[0],v0[1]+v1[1]];
        },
        sub: function(v0,v1) {
            return [v0[0]-v1[0],v0[1]-v1[1]];
        },
        mul: function(x,v) {
            return [x*v[0],x*v[1]];
        },
        copy: function(a) {
            var b = [];
            for (var i = 0; i < a.length; i++) {
                b[i] = a[i];
            }
            return b;
        },
        set: function(a,b) {
            for (var i = 0; i < b.length; i++) {
                a[i] = b[i];
            }
        }
    }

    // Basic partial differentiation
    function pdiff0(f,a,i,eps) {
        var ai = a[i];
        a[i] = ai-eps;
        var x0 = f(a);
        a[i] = ai+eps;
        var x1 = f(a);
        a[i] = ai;
        return (x1-x0)/(2*eps);
    }

    // Richardson extrapolation
    function pdiff1(f,a,i,eps)
    {
        var d1 = pdiff0(f,a,i,eps);
        var d2 = pdiff0(f,a,i,2*eps);
        return d1+(d1-d2)/3;
    }

    // Time derivatives
    function derivs0(f,p,v,dp,dv,t,eps)
    {
        var p1 = Vector.copy(p), v1 = Vector.copy(v);
        var p2 = Vector.copy(p), v2 = Vector.copy(v);
        f(p1,v1,t,eps);
        f(p2,v2,t,-eps);
        dp[0] = (p1[0]-p2[0])/(2*eps);
        dp[1] = (p1[1]-p2[1])/(2*eps);
        dv[0] = (v1[0]-v2[0])/(2*eps);
        dv[1] = (v1[1]-v2[1])/(2*eps);
    }

    function derivs1(f,p,v,dp,dv,t,eps)
    {
        var p1 = [], p2 = [], v1 = [], v2 = [];
        derivs0(f,p,v,p1,v1,t,eps/2);
        derivs0(f,p,v,p2,v2,t,eps);
        dp[0] = p1[0] + (p1[0]-p2[0])/3;
        dp[1] = p1[1] + (p1[1]-p2[1])/3;
        dv[0] = v1[0] + (v1[0]-v2[0])/3;
        dv[1] = v1[1] + (v1[1]-v2[1])/3;
    }

    // Choose your functions here.
    var derivs = derivs0; // Basic derivs
    var pdiff = pdiff0;   // Basic partials

    // d(df/dxj)/dxi - components of Hessian (Jacobian for PDs).
    function hessian(f,i,j,a,eps) {
        return pdiff(function(a) { return pdiff(f,a,j,eps); },
                     a,i,eps);
    }

    // Vector of partial derivs
    function pdiffs(f,a,eps)
    {
        var r = [];
        for (var i = 0; i < a.length; i++) {
            r[i] = pdiff(f,a,i,eps);
        }
        return r;
    }

    
    function single(f,p,v,t,dt,eps)
    {
        var dp = []
        var dv = []
        derivs(f,p,v,dp,dv,t,eps);
        p[0] += dt*dp[0];
        p[1] += dt*dp[1];
        v[0] += dt*dv[0];
        v[1] += dt*dv[1];
    }

    function rk4(f,p,v,t,dt,eps)
    {
        var p1 = [], v1 = [], p2 = [], v2 = [];
        var p3 = [], v3 = [], p4 = [], v4 = [];
        derivs(f,p,v,p1,v1,t,eps);
        derivs(f, 
               Vector.add(p,Vector.mul(0.5*dt,p1)),
               Vector.add(v,Vector.mul(0.5*dt,v1)),
               p2,v2,t,eps);
        derivs(f,
               Vector.add(p,Vector.mul(0.5*dt,p2)),
               Vector.add(v,Vector.mul(0.5*dt,v2)),
               p3,v3,t,eps);
        derivs(f,
               Vector.add(p,Vector.mul(dt,p3)),
               Vector.add(v,Vector.mul(dt,v3)),
               p4,v4,t,eps);
        p[0] += dt*(p1[0] + 2*p2[0] + 2*p3[0] + p4[0])/6;
        p[1] += dt*(p1[1] + 2*p2[1] + 2*p3[1] + p4[1])/6;
        v[0] += dt*(v1[0] + 2*v2[0] + 2*v3[0] + v4[0])/6;
        v[1] += dt*(v1[1] + 2*v2[1] + 2*v3[1] + v4[1])/6;
    }

    var steps = 1;
    function rk4var(f,p,v,t,dt,eps)
    {
        while (true) {
            var p1 = Vector.copy(p),v1 = Vector.copy(v);
            var p2 = Vector.copy(p),v2 = Vector.copy(v);
            rk4(f,p1,v1,t,dt/steps,eps);
            rk4(f,p2,v2,t,dt/(2*steps),eps);
            var p3 = Vector.copy(p2),v3 = Vector.copy(v2);
            rk4(f,p3,v3,t,dt/(2*steps),eps);
            //Sum the squares of the errors
            var errorvec = Vector.sub(v3,v1);
            var error = Vector.dot(errorvec,errorvec);
            if (error > 1e-7 && steps < 100) {
                steps *= 2;
                //console.log("Steps now " + steps);
                //alert("Steps now " + steps);
            } else if (error < 1e-6 && steps > 1) {
                steps /= 2;
                //alert("Steps now " + steps);
                //console.log("Steps now " + steps);
                Vector.set(p,p2); Vector.set(v,v2);
                break;
            } else {
                Vector.set(p,p3); Vector.set(v,v3);
                break;
            }
        }
        for (var i = 1; i < steps; i++) {
            rk4(f,p,v,t,dt/steps,eps);
        }
    }

    Lagrange.evaluate = function(L,p,v,t,eps) 
    {
        var pp = pdiffs(function(p) { return L(p,v,t) },
                        p,eps);
        var pv = pdiffs(function(v) { return L(p,v,t) },
                        v,eps);
        p[0] += v[0]*eps;
        p[1] += v[1]*eps;
        // Compute new values of pv (apply Euler-Lagrange)
        pv[0] += pp[0]*eps;
        pv[1] += pp[1]*eps;
        // pv is PD of L at t+eps
        // Now find new v' st pdiff(Lp,v') = pv
        // Initial guess is v, use Jacobian to refine guess
        // F(v') = F(v+dv) = F(v) + (JF)(dv), here F is the partial
        // derivatives of L at v (so it's really a Hessian).
        // Then v' = v + dv where dv = inv(JF)(F(v')-F(v))
        // Doing this step once seems to be enough
        // I'm surprised it works at all really.
        for (var k = 0; k < 1; k++) {
            var Lp = function(v) { return L(p,v,t) };
            var pv2 = pdiffs(Lp,v,eps);
            var t = Vector.sub(pv,pv2);
            var a = hessian(Lp,0,0,v,eps);
            var b = hessian(Lp,0,1,v,eps);
            var c = b; // Symmetry
            var d = hessian(Lp,1,1,v,eps);
            var det = a*d - b*c; // If det is too small, we'll just diverge here
            //process.stderr.write(a + " " + b + " " + c + " " + d + "\n");
            v[0] += (d*t[0] - b*t[1])/det;
            v[1] += (-c*t[0] + a*t[1])/det;
        }
    }

    Lagrange.step = rk4var;    // 4th order Runge Kutta, variable step size
    //var step = rk4;       // 4th order Runge Kutta
    //var step = single;    // This one isn't too good.

    
    // Some Lagrangians
    Lagrange.OrbitPolar = function(G)
    {
        this.L = function(p,v,t) {
            var m = 1;
            var r = p[0];
            var rv = v[0];
            var thetav = v[1];
            var ke = m*0.5*(rv*rv + r*thetav*r*thetav);
            var pe = -m*G/r;
            var l = ke-pe;
            return l;
        }
        this.coords = function(p,v) {
            return [p[0]*Math.sin(p[1]), p[0]*Math.cos(p[1]),0,0,0,0];
        }
    }

    Lagrange.OrbitRectangular = function(G)
    {
        this.L = function(p,v,t) {
            var m = 1;
            var r = Vector.length(p);
            var ke = m*0.5*(Vector.dot(v,v));
            var pe = -m*G/r;
            var l = ke-pe;
            return l;
        }
        this.coords = function(p,v,t) {
            return [p[1],p[0],0,0,0,0];
        }
    }

    Lagrange.DoubleOrbit = function(G)
    {
        this.L = function(p,v,t) {
            var m1 = 1, m2 = 1;
            var r1 = [p[0],p[1]];
            var r2 = [p[0]+10,p[1]];
            var pe = -(m1*G/Vector.length(r1) + m2*G/Vector.length(r2));
            var ke = 0.5*(Vector.dot(v,v));
            var l = ke-pe;
            return l;
        }
        this.coords = function(p,v,t) {
            return [p[1],p[0],0,0,0,0];
        }
    }

    Lagrange.Hookean = function(H,m) {
        this.L = function(p,v,t) {
            var pe = H*Vector.length(p);
            var ke = 0.5*(Vector.dot(v,v));
            var l = ke - pe;
            return l;
        }
        this.coords = function(p,v,t) {
            return [p[1],p[0],0];
        }
    }

    Lagrange.Hookean2 = function(H,m) {
        this.L = function(p,v,t) {
            var p1 = [p[0],p[1]];
            var p2 = [p[0]+0.5,p[1]];
            var pe1 = H*Vector.length(p1);
            var pe2 = H*Vector.length(p2);
            var ke = 0.5*(Vector.dot(v,v));
            var l = ke - pe1 -pe2;
            return l;
        }
        this.coords = function(p,v,t) {
            return [p[1],p[0],0];
        }
    }

    Lagrange.SphericalPendulum = function(G,r)
    {
        this.L = function(p,v,t) {
            var m = 1;
            var theta = p[0];
            var phi = p[1];
            var thetav = v[0];
            var phiv = v[1];
            var ke = 0.5*m*r*r*(Math.pow(thetav,2) + Math.pow(Math.sin(theta)*phiv,2));
            var pe = -m*G*r*Math.cos(theta);
            var l = ke-pe;
            return l;
        }
        this.coords = function(p,v,t) {
            var theta = p[0];
            var phi = p[1];
            var z = r* Math.sin(phi)*Math.sin(theta);
            var y = -r*Math.cos(theta);
            var x = r*Math.cos(phi)*Math.sin(theta);
            return [x,y,z,theta,phi,0];
        }
    }

    Lagrange.DoublePendulum = function(G,m0,m1,r0,r1) {
        this.L = function(p,v,t) {
            var y0,y1,v0,v1,ke0,ke1,pe0,pe1;
            y0 = r0*Math.sin(p[0]);
            y1 = y0 + r1*Math.sin(p[1]);
            v0 = [-r0*Math.sin(p[0])*v[0], r0*Math.cos(p[0])*v[0]];
            v1 = [-r1*Math.sin(p[1])*v[1], r1*Math.cos(p[1])*v[1]];
            ke0 = 0.5*m0*(Vector.dot(v0,v0));
            ke1 = 0.5*m1*(Vector.dot(Vector.add(v0,v1),
                                     Vector.add(v0,v1)));
            pe0 = G*m0*y0;
            pe1 = G*m1*y1;
            return ke0+ke1-pe0-pe1;
        }
        this.coords = function(p,v,t) {
            var x0,y0,x1,y1;
            x0 = r0*Math.cos(p[0]);
            y0 = r0*Math.sin(p[0]);
            x1 = x0 + r1*Math.cos(p[1]);
            y1 = y0 + r1*Math.sin(p[1]);
            return [x1,y1,0,x0,y0,0];
        }
    }

    Lagrange.run = function(system,p,v,dt,eps) {
        var f = function(p,v,t,eps) { evaluate(system.L,p,v,t,eps); }
        console.log(system.coords(p,v).join(" "));
        for (var i = 0; i < 10000; i++) {
            var t = i*dt;
            step(f,p,v,t,dt,eps);
            console.log(system.coords(p,v).join(" "));
        }
    }
}())

//run(new OrbitPolar(0.0002),[1,0],[0.0,0.01],1,1e-4);
//run(new OrbitRectangular(0.0002),[1,0],[0.0,0.01],1,1e-4);
//run(new DoubleOrbit(0.0002),[1,0],[0.0,0.01],1,1e-4);
//run(new Hookean(0.001,1),[1,0],[0.0,0.02],1,1e-4);
//run(new Hookean2(0.001,1),[1,0],[0.0,0.02],1,1e-4);
//run(new SphericalPendulum(0.0002,1), [1,0],[0.01,0.05],1,1e-4);
//run(new DoublePendulum(0.0002,1,1,0.5,0.5),[1,0],[0.0,0.01],1,1e-4);
