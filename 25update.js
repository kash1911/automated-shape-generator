/*******************************************************
 * p5.js + Dubins Path + 出界检测 + TSP + 5~8点限制
 *
 * 1) 在 draw() 中做随机网格Dubins绘制：
 *    - 若路径越界，则不画。
 *    - 同时把点收集到 allPoints。
 * 2) TSP (限5~8点) 也用出界检测：
 *    - 若越界 => cost=∞
 *    - 否则累加离散距离
 * 3) 最后绘制 Dubins 闭合环(首尾也Dubins)
 *******************************************************/

let pwidth = 400;
let pheight = 400;

let num = 5;
let current = 0;
let order = 0;
let allPoints = []; // 存所有“实际使用”的点

function setup() {
  createCanvas(800, 800, SVG);
  colorMode(RGB, 100);
  // drawingContext.shadowColor  = color(0,0,0,15);
  // drawingContext.shadowBlur = width / 10;
}

function draw() {
  background(255);
  strokeWeight(8);

  let origin_x, origin_y, origin_yaw;
  let start_x, start_y, start_yaw;
  let end_x, end_y, end_yaw;
  let curvature;

  // 下面是原先用于随机的数组
  let nu = [1, -1, 0];
  let n = [1, -1];
  let ny = [0, -1];
  let ang = [0, HALF_PI, PI, 3 * HALF_PI, TWO_PI];
  let ranu = random(nu);
  let ran = random(n);
  let rapy = random(ny);
  let rangle = random(ang);

  let tilesX = 4;
  let tilesY = tilesX;
  let gap = 15;

  let tileW = pwidth / tilesX;
  let tileH = pheight / tilesY;
  curvature = (tileH - gap) / 2;

  /************************************
   * 1) 原始随机绘制
   ************************************/
  //beginShape();

  for (let x = 0; x < tilesX; x++) {
    for (let y = 0; y < tilesY; y++) {
      let posX = tileW * x - gap;
      let posY = tileH * y - gap;
      let randomSquare = floor(random(2));

      if (randomSquare == 0) {
        order++;
        if (order === 1) {
          // 第一次
          start_x = posX + (tileW / 2) + gap + ranu * (tileW - gap) / 2;
          if (ranu === 0) {
            if (rapy === -1) {
              start_y = tileH * y - rapy * gap / 2;
            } else {
              start_y = tileH * y + tileH - gap / 2;
            }
          } else {
            start_y = posY + (tileH / 2) + gap;
          }
          start_yaw = rangle;

          origin_x = start_x;
          origin_y = start_y;
          origin_yaw = start_yaw;

          allPoints.push({ x: start_x, y: start_y, yaw: start_yaw });
        }

        end_x = posX + (tileW / 2) + gap + ranu * (tileW - gap) / 2;
        if (ranu === 0) {
          if (rapy === -1) {
            end_y = tileH * y - rapy * gap / 2;
          } else {
            end_y = tileH * y + tileH - gap / 2;
          }
        } else {
          end_y = posY + (tileH / 2) + gap;
        }
        end_yaw = start_yaw * ran;

        // 改动：用 buildDiscreteDubins() 出界检测
        let seg = buildDiscreteDubins(
          { x: start_x, y: start_y, yaw: start_yaw },
          { x: end_x,   y: end_y,   yaw: end_yaw },
          curvature
        );
        if (seg) {
          // 没有越界 => 绘制
          for (let i = 0; i < seg.length; i++) {
            vertex(seg[i].x, seg[i].y);
          }
          // 并收集点
          allPoints.push({ x: end_x, y: end_y, yaw: end_yaw });

          // 更新
          start_x = end_x;
          start_y = end_y;
          start_yaw = end_yaw;
        } else {
          // 出界 => 不画 & 不更新
          console.log("随机路径越界，跳过");
        }
      }
    }
  }
  // 回到起点(用同样方式检测越界)
  let backSeg = buildDiscreteDubins(
    { x: start_x, y: start_y, yaw: start_yaw },
    { x: origin_x, y: origin_y, yaw: origin_yaw },
    curvature
  );
  if (backSeg) {
    for (let i = 0; i < backSeg.length; i++) {
      vertex(backSeg[i].x, backSeg[i].y);
    }
    // 收集最终点
    allPoints.push({ x: origin_x, y: origin_y, yaw: origin_yaw });
  } else {
    console.log("最后一段回到起点越界，跳过");
  }

  //endShape();

  /************************************
   * 2) TSP 优化 + 5~8 点限制
   ************************************/
  if (allPoints.length < 5) {
    console.log("收集到点数=", allPoints.length, "不足5个, 不做TSP");
  } else {
    let desiredCount = floor(random(6, 18)); // 5..8
    let finalPoints = [];
    if (allPoints.length <= desiredCount) {
      finalPoints = allPoints;
    } else {
      finalPoints = randomSubset(allPoints, desiredCount);
    }

    // 构建代价矩阵(出界 => cost=∞)
    let costMat = buildCostMatrix(finalPoints, curvature);
    let res = tspHeldKarp(costMat);
    let bestOrder = res.bestOrder;
    let bestCost = res.bestCost;

    if (!bestOrder) {
      console.log("TSP无可行解");
    } else {
      console.log("TSP: bestOrder=", bestOrder, "bestCost=", bestCost);

      // 画 TSP 路径(绿色)
      stroke(0, 150, 0);
      strokeWeight(3);
      noFill();

      beginShape();

      // (a) 中间段 i-> i+1
      for (let i = 0; i < bestOrder.length - 1; i++) {
        let i1 = bestOrder[i];
        let i2 = bestOrder[i + 1];
        let seg2 = buildDiscreteDubins(finalPoints[i1], finalPoints[i2], curvature);
        if (seg2) {
          for (let pt of seg2) {
            vertex(pt.x, pt.y);
          }
        }
      }
      // (b) 最后一段 last->first
      let lastIndex = bestOrder[bestOrder.length - 1];
      let firstIndex = bestOrder[0];
      let seg3 = buildDiscreteDubins(finalPoints[lastIndex], finalPoints[firstIndex], curvature);
      if (seg3) {
        for (let pt of seg3) {
          vertex(pt.x, pt.y);
        }
      }

      endShape();
      save("mySVG.svg");
    print("saved svg");
    }
  }

  noLoop();
  //if (mouseIsPressed) {
    //save("mySVG.svg");
    //print("saved svg");
  //}
}

/*****************************************
 * 出界检测：buildDiscreteDubins()
 * 
 * 若任何离散点不在[0..pwidth]×[0..pheight],
 * 则返回null。否则返回离散点数组
 *****************************************/
function buildDiscreteDubins(p1, p2, c) {
  let dub = dubins_path_planning(
    p1.x, p1.y, p1.yaw,
    p2.x, p2.y, p2.yaw,
    c
  );
  if (!dub) return null;
  let px = dub[0];
  let py = dub[1];
  if (px.length < 2) return null;

  let arr = [];
  for (let i = 0; i < px.length; i++) {
    let xx = px[i];
    let yy = py[i];
    // 出界检测
    if (xx < 0 || xx > pwidth || yy < 0 || yy > pheight) {
      return null; // 直接返回null => 放弃这条路径
    }
    arr.push({ x: xx, y: yy });
  }
  return arr;
}

/*****************************************
 * 其余Dubins函数(与之前一致)
 *****************************************/
function mod2pi(theta) {
  return theta - 2 * PI * floor(theta / (2 * PI));
}
function pi_2_pi(angle) {
  while (angle >= PI) angle -= 2 * PI;
  while (angle <= -PI) angle += 2 * PI;
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
  let dx = ex;
  let dy = ey;
  let e = sqrt(dx ** 2 + dy ** 2);
  let d = e / c;

  let theta = mod2pi(atan2(dy, dx));
  let alpha = mod2pi(-theta);
  let beta = mod2pi(eyaw - theta);

  let planners = [LSL, RSR, LSR, RSL, RLR, LRL];
  let bcost = Infinity;
  let bt, bp, bq, bmode;
  for (planner of planners) {
    let arr = planner(alpha, beta, d);
    let t = arr[0];
    let p = arr[1];
    let q = arr[2];
    let mode = arr[3];
    if (t == null) {
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
    let d;
    if (m == "S") {
      d = 1.0 / c;
    } else {
      d = radians(1.0 / 2);
    }
    while (pd < abs(l - d)) {
      let pxv = px[px.length - 1] + d * c * cos(pyaw[pyaw.length - 1]);
      let pyv = py[py.length - 1] + d * c * sin(pyaw[pyaw.length - 1]);
      px.push(pxv);
      py.push(pyv);
      if (m == "L") {
        pyaw.push(pyaw[pyaw.length - 1] + d);
      } else if (m == "S") {
        pyaw.push(pyaw[pyaw.length - 1]);
      } else if (m == "R") {
        pyaw.push(pyaw[pyaw.length - 1] - d);
      }
      pd += d;
    }
    d = l - pd;
    px.push(px[px.length - 1] + d * c * cos(pyaw[pyaw.length - 1]));
    py.push(py[py.length - 1] + d * c * sin(pyaw[pyaw.length - 1]));

    if (m == "L") {
      pyaw.push(pyaw[pyaw.length - 1] + d);
    } else if (m == "S") {
      pyaw.push(pyaw[pyaw.length - 1]);
    } else if (m == "R") {
      pyaw.push(pyaw[pyaw.length - 1] - d);
    }
    pd += d;
  }
  return [px, py, pyaw];
}

function dubins_path_planning(sx, sy, syaw, ex, ey, eyaw, c) {
  ex = ex - sx;
  ey = ey - sy;

  let lex = cos(syaw) * ex + sin(syaw) * ey;
  let ley = -sin(syaw) * ex + cos(syaw) * ey;
  let leyaw = eyaw - syaw;
  let arr = dubins_path_planning_from_origin(lex, ley, leyaw, c);

  let lpx = arr[0];
  let lpy = arr[1];
  let lpyaw = arr[2];
  let mode = arr[3];
  let clen = arr[4];

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
    let iyaw = lpyaw[i];
    pyaw.push(pi_2_pi(iyaw + syaw));
  }

  return [px, py, pyaw, mode, clen];
}

/*****************************************
 * TSP + 辅助
 *****************************************/
function buildCostMatrix(pts, c) {
  let N = pts.length;
  let mat = new Array(N).fill(null).map(()=> new Array(N).fill(Infinity));
  for (let i=0; i<N; i++){
    for (let j=0; j<N; j++){
      if(i===j){
        mat[i][j]= Infinity;
      } else {
        let seg = buildDiscreteDubins(pts[i], pts[j], c);
        if(!seg){
          mat[i][j]= Infinity;
        } else {
          let distSum=0;
          for(let k=1; k<seg.length; k++){
            distSum+= dist(seg[k].x, seg[k].y, seg[k-1].x, seg[k-1].y);
          }
          mat[i][j]= distSum;
        }
      }
    }
  }
  return mat;
}

function tspHeldKarp(costMat){
  let N= costMat.length;
  if(N<2) return {bestCost:Infinity, bestOrder:null};

  let dp= new Array(1<<N).fill(null).map(()=> new Array(N).fill(Infinity));
  let parent= new Array(1<<N).fill(null).map(()=> new Array(N).fill(-1));

  dp[1][0]=0; // 从0出发

  for(let mask=1; mask<(1<<N); mask++){
    for(let i=0; i<N; i++){
      if((mask & (1<<i))===0) continue;
      let baseCost= dp[mask][i];
      if(baseCost===Infinity) continue;
      for(let j=0; j<N; j++){
        if(j===i) continue;
        if((mask&(1<<j))!==0) continue;
        let nm= mask|(1<<j);
        let newCost= baseCost+ costMat[i][j];
        if(newCost< dp[nm][j]){
          dp[nm][j]= newCost;
          parent[nm][j]= i;
        }
      }
    }
  }

  let fullMask= (1<<N)-1;
  let bestCost= Infinity;
  let bestLast= -1;
  for(let i=1; i<N; i++){
    let cst= dp[fullMask][i]+ costMat[i][0];
    if(cst< bestCost){
      bestCost= cst; bestLast= i;
    }
  }
  if(bestCost===Infinity) return {bestCost:Infinity, bestOrder:null};

  // 回溯
  let route=[];
  let cur= bestLast;
  let cm= fullMask;
  while(cur!==-1){
    route.push(cur);
    let p= parent[cm][cur];
    cm= cm^(1<<cur);
    cur= p;
  }
  route.push(0);
  route.reverse();

  return {bestCost, bestOrder: route};
}

// 若抽取点>8需要随机挑
function randomSubset(arr, n){
  if(n>= arr.length) return arr;
  let copy= [...arr];
  shuffle(copy, true); // p5.js内置 shuffle
  return copy.slice(0, n);
}
