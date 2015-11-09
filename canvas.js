"use strict"

function simulate(canvas,system)
{
    //var args = location.search;
    //alert("'" + args + "'");
    var p = [1,0], v = [0,1], eps = 1e-4;

    // Time is ms so
    function run1(system,t,dt) {
        var maxstep = 5
        while (dt > 0) {
            var f = function(p,v,t,eps) { evaluate(system.L,p,v,t,eps); }
            var delta = Math.min(maxstep,dt)
            step(f,p,v,t,delta/1000,eps);
            dt -= delta
        }
        return system.coords(p,v);
    }

    var ctx = canvas.getContext("2d");
    var N = 1024;
    var M = 64;

    var lasttime = 0;

    var xs = [];
    var ys = [];
    var start = 0;
    var psize = 0;
    var maxpsize = 150;

    function redisplay2(timestamp) {
        console.log(timestamp);
        if (lasttime == 0) {
            lasttime = timestamp;
        }
        var coords = run1(system,timestamp,timestamp-lasttime)
        lasttime = timestamp;
        ctx.clearRect(0,0,800,800);
        var ox = 400, oy = 300
        var l1 = 250, l2 = 200;
        var s1 = 5000;
        var s2 = 1234;
        var x1 = ox + l1*coords[3];
        var y1 = oy - l1*coords[4];
        var x2 = ox + l1*coords[0];
        var y2 = oy - l1*coords[1];
        //var x1 = ox+l1*Math.sin(timestamp/s1);
        //var y1 = oy+l1*Math.cos(timestamp/s1);
        //var x2 = x1 + l2*Math.sin(timestamp/s2);
        //var y2 = y1 + l2*Math.cos(timestamp/s2);
        xs[(start+psize)%maxpsize] = x2;
        ys[(start+psize)%maxpsize] = y2;
        if (psize < maxpsize-1) psize++;
        else start++;

        ctx.lineStyle = "#00FF00";
        //ctx.fillStyle = "#00FF00";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(xs[start],ys[start]);
        for (var i = 1; i < psize; i++) {
            ctx.lineTo(xs[(start+i)%maxpsize],ys[(start+i)%maxpsize]);
        }
        ctx.stroke();

        ctx.lineWidth = 5;
        ctx.fillStyle = "#FF0000";
        ctx.beginPath();
        ctx.moveTo(ox,oy);
        ctx.lineTo(x1,y1);
        ctx.lineTo(x2,y2);
        ctx.stroke();
        ctx.beginPath();
        ctx.arc(ox,oy,10,0,2*Math.PI);
        ctx.fill();
        ctx.beginPath();
        ctx.arc(x1,y1,10,0,2*Math.PI);
        ctx.fill();
        ctx.beginPath();
        ctx.arc(x2,y2,10,0,2*Math.PI);
        ctx.fill();
        window.requestAnimationFrame(redisplay2);
    }

    window.requestAnimationFrame(redisplay2);
}

