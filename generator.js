let num =5;
let current = 0;
let order = 0;



function setup() {
  createCanvas(800, 800,SVG);
  colorMode(RGB, 100);
  
  // drawingContext.shadowColor  = color(0,0,0,15);
  // drawingContext.shadowBlur = width / 10;
}

function draw() {
  background(255);
  
  strokeWeight(8);
  // noStroke();
  // fill(0);
  
  let origin_x, origin_y, origin_yaw;
  
    
    // let fc = color(0);
  let start_x;
  let start_y;
  let start_yaw;
  let end_x;
  let end_y;
  let end_yaw;
  let curvature;
  let arrs;
   let pwidth = 400;
  let pheight = 400;
  
  let nu = [1,-1,0];
    let n = [1,-1];
   let ny = [0,-1];
    let ang = [0,HALF_PI,PI,3*HALF_PI,TWO_PI];
    let ranu = random(nu);
    let ran = random(n);
  let rapy = random(ny);
    let rangle = random(ang);
  
  

  let tilesX = 4;
  let tilesY = tilesX;
  let gap = 0;
  
  let tileW = pwidth / tilesX;
  let tileH = pheight / tilesY;
  curvature = (tileH - gap)/2;
  
  
  
  
  beginShape();
  push();
  
  translate(pwidth/2,pheight/2);
  for (let x = 0; x < tilesX; x++) {
    for (let y = 0; y < tilesY; y++) {
    
    let posX = tileW * x - gap;
    let posY = tileH * y - gap;
    let randomSquare = floor(random(2));
    
    // push();
    
    // translate(posX +(tileW/2), posY+ (tileH/2));
    
    if (randomSquare == 0) {
      // fill(fc);
      order ++;
      
      
      
      if (order == 1) {
        
          
        
     start_x = posX +(tileW/2) +gap + ranu * (tileW - gap)/2 ;
         if(ranu == 0){
                if(rapy == -1 ){
                 start_y =  tileH * y - rapy * gap/2;
                }else{
                    
                  start_y =  tileH * y + tileH - gap/2;
                  
                }
                 }else{
                 
              start_y = posY+ (tileH/2) + gap;
            }
      
      start_yaw = rangle;
      origin_x = start_x;
      origin_y = start_y;
      origin_yaw = start_yaw;
         
    }
   end_x = posX +(tileW/2) +gap + ranu * (tileW - gap)/2 ;
         if(ranu == 0){
                if(rapy == -1 ){
                 end_y =  tileH * y - rapy * gap/2;
                }else{
                    
                  end_y =  tileH * y + tileH - gap/2;
                  
                }
                 }else{
                 
              end_y = posY+ (tileH/2) + gap;
            }
    end_yaw = start_yaw * ran;
      
      
    let arrs = dubins_path_planning(
      start_x,
      start_y,
      start_yaw,
      end_x,
      end_y,
      end_yaw,
      curvature
    );
    let px = arrs[0];
    let py = arrs[1];

      
      
    for (let i = 0; i < px.length; i++) {
      let x = px[i];
      let y = py[i];
          
      vertex(x, y);
    }
    start_x = end_x;
    start_y = end_y;
    start_yaw = end_yaw;

  
      
    // } else {
    //   noFill();
    }
    
    
      
      
    
    // pop();
    
    }
  }
  
  
  

  
  

      
      
    
  

  arrs = dubins_path_planning(
    end_x,
    end_y,
    end_yaw,
    origin_x,
    origin_y,
    origin_yaw,
    curvature
  );
  let px = arrs[0];
  let py = arrs[1];

  for (let i = 0; i < px.length; i++) {
      
      push();
  
      translate(width/2,height/2);
  
      
    let x = px[i];
    let y = py[i];
    vertex(x, y);
      
     pop(); 
  }
  // ellipse(origin_x, origin_y, 5, 5);
  endShape();

  noLoop();
  if(mouseIsPressed === true) {
     save("mySVG.svg"); // give file name
  print("saved svg");
  }
     
  // frameRate(1);
}

//0〜TWO_PIで値を返す
function mod2pi(theta) {
  return theta - 2.0 * PI * floor(theta / 2.0 / PI);
}

function pi_2_pi(angle) {
  while (angle >= PI) {
    angle = angle - TWO_PI;
  }
  while (angle <= -PI) {
    angle = angle + TWO_PI;
  }
  return angle;
}

function LSL(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);
  let tmp0 = d + sa - sb;

  let mode = ["L", "S", "L"];
  let p_squared = 2 + d * d - 2 * c_ab + 2 * d * (sa - sb);
  if (p_squared < 0) {
    return [null, null, null, mode];
  }
  let tmp1 = atan2(cb - ca, tmp0);
  let t = mod2pi(-alpha + tmp1);
  let p = sqrt(p_squared);
  let q = mod2pi(beta - tmp1);

  return [t, p, q, mode];
}

function RSR(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);

  let tmp0 = d - sa + sb;
  let mode = ["R", "S", "R"];
  let p_squared = 2 + d * d - 2 * c_ab + 2 * d * (sb - sa);
  if (p_squared < 0) {
    return [null, null, null, mode];
  }
  let tmp1 = atan2(ca - cb, tmp0);
  let t = mod2pi(alpha - tmp1);
  let p = sqrt(p_squared);
  let q = mod2pi(-beta + tmp1);

  return [t, p, q, mode];
}

function LSR(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);

  let p_squared = -2 + d * d + 2 * c_ab + 2 * d * (sa + sb);
  let mode = ["L", "S", "R"];
  if (p_squared < 0) {
    return [null, null, null, mode];
  }
  let p = sqrt(p_squared);
  let tmp2 = atan2(-ca - cb, d + sa + sb) - atan2(-2.0, p);
  let t = mod2pi(-alpha + tmp2);
  let q = mod2pi(-mod2pi(beta) + tmp2);

  return [t, p, q, mode];
}

function RSL(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);

  let p_squared = d * d - 2 + 2 * c_ab - 2 * d * (sa + sb);
  let mode = ["R", "S", "L"];
  if (p_squared < 0) {
    return [null, null, null, mode];
  }
  let p = sqrt(p_squared);
  let tmp2 = atan2(ca + cb, d - sa - sb) - atan2(2.0, p);
  let t = mod2pi(alpha - tmp2);
  let q = mod2pi(beta - tmp2);
  return [t, p, q, mode];
}

function RLR(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);

  let mode = ["R", "L", "R"];
  let tmp_rlr = (6.0 - d * d + 2.0 * c_ab + 2.0 * d * (sa - sb)) / 8.0;
  if (abs(tmp_rlr) > 1.0) {
    return [null, null, null, mode];
  }

  let p = mod2pi(2 * PI - acos(tmp_rlr));
  let t = mod2pi(alpha - atan2(ca - cb, d - sa + sb) + mod2pi(p / 2.0));
  let q = mod2pi(alpha - beta - t + mod2pi(p));
  return [t, p, q, mode];
}

function LRL(alpha, beta, d) {
  let sa = sin(alpha);
  let sb = sin(beta);
  let ca = cos(alpha);
  let cb = cos(beta);
  let c_ab = cos(alpha - beta);

  let mode = ["L", "R", "L"];
  let tmp_lrl = (6 - d * d + 2 * c_ab + 2 * d * (-sa + sb)) / 8;
  if (abs(tmp_lrl) > 1.0) {
    return [null, null, null, mode];
  }
  let p = mod2pi(2 * PI - acos(tmp_lrl));
  let t = mod2pi(-alpha - atan2(ca - cb, d + sa - sb) + p / 2);
  let q = mod2pi(mod2pi(beta) - alpha - t + mod2pi(p));

  return [t, p, q, mode];
}

function dubins_path_planning_from_origin(ex, ey, eyaw, c) {
  //# nomalize
  let dx = ex;
  let dy = ey;
  let e = sqrt(dx ** 2 + dy ** 2);
  let d = e / c;
  // print(dx, dy, e, d)

  let theta = mod2pi(atan2(dy, dx));
  let alpha = mod2pi(-theta);
  let beta = mod2pi(eyaw - theta);

  // print(theta, alpha, beta, d)

  let planners = [LSL, RSR, LSR, RSL, RLR, LRL];
  let bcost = Infinity;
  let bt;
  let bp;
  let bq;
  let bmode;
  for (planner of planners) {
    // print(planner);
    let arr = planner(alpha, beta, d);
    // print(arr);
    let t = arr[0];
    let p = arr[1];
    let q = arr[2];
    let mode = arr[3];
    if (t == null) {
      // print(mode + " cannot generate path");
      continue;
    }

    let cost = abs(t) + abs(p) + abs(q);
    if (bcost > cost) {
      bt = t;
      bp = p;
      bq = q;
      bmode = mode;
      bcost = cost;
    }
  }

  let arr2 = generate_course([bt, bp, bq], bmode, c);

  let px = arr2[0];
  let py = arr2[1];
  let pyaw = arr2[2];
  // print(px, py, pyaw, bmode, bcost);

  return [px, py, pyaw, bmode, bcost];
}

function generate_course(length, mode, c) {
  let px = [0];
  let py = [0];
  let pyaw = [0];

  for (let i = 0; i < length.length; i++) {
    let m = mode[i];
    let l = length[i];
    let pd = 0.0;
    if (m == "S") {
      d = 1.0 / c;
    } else {
      //# turning couse
      d = radians(1.0 / 2);
    }
    while (pd < abs(l - d)) {
      let pxv = px[px.length - 1] + d * c * cos(pyaw[pyaw.length - 1]);
      let pyv = py[py.length - 1] + d * c * sin(pyaw[pyaw.length - 1]);
      px.push(pxv);
      py.push(pyv);
      if (m == "L") pyaw.push(pyaw[pyaw.length - 1] + d);
      else if (m == "S") pyaw.push(pyaw[pyaw.length - 1]);
      else if (m == "R") pyaw.push(pyaw[pyaw.length - 1] - d);
      pd += d;
    }
    d = l - pd;
    px.push(px[px.length - 1] + d * c * cos(pyaw[pyaw.length - 1]));
    py.push(py[py.length - 1] + d * c * sin(pyaw[pyaw.length - 1]));

    if (m == "L") pyaw.push(pyaw[pyaw.length - 1] + d);
    else if (m == "S") pyaw.push(pyaw[pyaw.length - 1]);
    else if (m == "R") pyaw.push(pyaw[pyaw.length - 1] - d);
    pd += d;
  }
  return [px, py, pyaw];
}

function dubins_path_planning(sx, sy, syaw, ex, ey, eyaw, c) {
  // """
  // Dubins path plannner
  // input:
  //     sx x position of start point [m]
  //     sy y position of start point [m]
  //     syaw yaw angle of start point [rad]
  //     ex x position of end point [m]
  //     ey y position of end point [m]
  //     eyaw yaw angle of end point [rad]
  //     c curvature [1/m]
  // output:
  //     px
  //     py
  //     pyaw
  //     mode
  // """

  ex = ex - sx;
  ey = ey - sy;

  let lex = cos(syaw) * ex + sin(syaw) * ey;
  let ley = -sin(syaw) * ex + cos(syaw) * ey;
  let leyaw = eyaw - syaw;
  let arr = dubins_path_planning_from_origin(lex, ley, leyaw, c);
  // print(arr);
  let lpx, lpy, lpyaw, mode, clen;
  lpx = arr[0];
  lpy = arr[1];
  lpyaw = arr[2];
  mode = arr[3];
  clen = arr[4];

  let px = [];
  let py = [];
  let pyaw = [];

  for (let i = 0; i < lpx.length; i++) {
    let x = lpx[i];
    let y = lpy[i];
    px.push(cos(-syaw) * x + sin(-syaw) * y + sx);
    py.push(-sin(-syaw) * x + cos(-syaw) * y + sy);
  }
  for (let i = 0; i < lpyaw.length; i++) {
    iyaw = lpyaw[i];
    pyaw.push(pi_2_pi(iyaw + syaw));
  }
  //     #  print(syaw)
  //     #  pyaw = lpyaw

  //     #  plt.plot(pyaw, "-r")
  //     #  plt.plot(lpyaw, "-b")
  //     #  plt.plot(eyaw, "*r")
  //     #  plt.plot(syaw, "*b")
  //     #  plt.show()

  return [px, py, pyaw, mode, clen];
}
