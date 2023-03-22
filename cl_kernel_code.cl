//Preprocessor things for compilation of tnp
#ifndef XLOW
  #define XLOW 0.f
#endif
#ifndef YLOW
  #define YLOW 0.f
#endif
#ifndef ZLOW
  #define ZLOW 0.f
#endif
#ifndef XHIGH
  #define XHIGH 0.f
#endif
#ifndef YHIGH
  #define YHIGH 0.f
#endif
#ifndef ZHIGH
  #define ZHIGH 0.f
#endif
#ifndef DX
  #define DX 0
#endif
#ifndef DY
  #define DY 0
#endif
#ifndef DZ
  #define DZ 0
#endif
#ifndef NX
  #define NX 0
#endif
#ifndef NY
  #define NY 0
#endif
#ifndef NZ
  #define NZ 0
#endif

void kernel vector_cross_mul(global float *A0, global const float *B0,
                             global const float *C0, global float *A1,
                             global const float *B1, global const float *C1,
                             global float *A2, global const float *B2,
                             global const float *C2) {
  int i = get_global_id(0); // Get index of the current element to be processed
  A0[i] = B1[i] * C2[i] - B2[i] * C1[i]; // Do the operation
  A1[i] = B2[i] * C0[i] - B0[i] * C2[i];
  A2[i] = B0[i] * C1[i] - B1[i] * C0[i];
}

void kernel vector_mul(global float *A, global const float *B,
                       global const float *C) {
  int i = get_global_id(0); // Get index of the current element to be processed
  A[i] = B[i] * C[i];       // Do the operation
}

void kernel vector_muls_addv(global float *A, global const float *B,
                             global const float *C) {
  float Bb = B[0];
  int i = get_global_id(0); // Get index of current element processed
  A[i] = Bb * A[i] + C[i];  // Do the operation
}

void kernel vector_muls(global float *A, global const float *B) {
  float Bb = B[0];
  int i = get_global_id(0); // Get index of current element processed
  A[i] = Bb * A[i];         // Do the operation
}

void kernel vector_mul_complex(global float2 *A, global float2 *B, global float2 *C) {
  int i = get_global_id(0); // Get index of the current element to be processed
  float2 b = B[i], c = C[i];
  A[i] = (float2)(b.s0 * c.s0 - b.s1 * c.s1, b.s0 * c.s1 + b.s1 * c.s0);
}

void kernel
tnp_k_optimized(global const float8 *a1, global const float8 *a2,                   // E, B coeff
      global float *x0, global float *y0, global float *z0, // initial pos
      global float *x1, global float *y1, global float *z1,
      const float Bcoeff, const float Ecoeff, // Bcoeff, Ecoeff
      const unsigned int n, const unsigned int ncalc // n, ncalc
      ) {
  uint id = get_global_id(0);
  uint prev_idx = UINT_MAX;
  float xprev = x0[id], yprev = y0[id], zprev = z0[id],
        x = x1[id], y = y1[id], z = z1[id];
  float8 temp, pos;
  float8 store0, store1, store2, store3, store4, store5;
  for (int t = 0; t < ncalc; t++) {
    float xy = x * y, xz = x * z, yz = y * z, xyz = x * yz;
    uint idx = ((uint)((z - ZLOW) / DZ) * NZ + (uint)((y - YLOW) / DY)) * NY + (uint)((x - XLOW) / DX);
    // round down the cells - this is intentional
    idx *= 3;
    pos = (float8)(1.f, x, y, z, xy, xz, yz, xyz);

    // Is there no better way to do this? Why does float8 not have dot()?
    if (prev_idx != idx) {
      store0 = a1[idx]; store1 = a1[idx+1]; store2 = a1[idx+2]; store3 = a2[idx]; store4 = a2[idx+1]; store5 = a2[idx+2];
      prev_idx = idx;
    }
    temp = store0 * pos;
    float xE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store1 * pos;
    float yE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store2 * pos;
    float zE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store3 * pos;
    float xP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store4 * pos;
    float yP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store5 * pos;
    float zP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;

    xP *= Bcoeff;
    yP *= Bcoeff;
    zP *= Bcoeff;
    xE *= Ecoeff;
    yE *= Ecoeff;
    zE *= Ecoeff;

    float xyP = xP * yP, yzP = yP * zP, xzP = xP * zP;
    float xxP = xP * xP, yyP = yP * yP, zzP = zP * zP;
    float r_det = 1.0 / (1.0 + xxP + yyP + zzP);

    float dx = 2.0 * x - xprev - zP * yprev + yP * zprev;
    float dy = 2.0 * y - yprev - xP * zprev + zP * xprev;
    float dz = 2.0 * z - zprev - yP * xprev + xP * yprev;
  
    xprev = x;
    yprev = y;
    zprev = z;

    x = clamp(xE + r_det * ((1.f + xxP) * dx + (zP + xyP) * dy + (xzP - yP) * dz), XLOW, XHIGH);
    y = clamp(yE + r_det * ((xyP - zP) * dx + (1.f + yyP) * dy + (yzP + xP) * dz), YLOW, YHIGH);
    z = clamp(zE + r_det * ((xzP + yP) * dx + (yzP - xP) * dy + (1.f + zzP) * dz), ZLOW, ZHIGH);
    //if(id == 0) printf("%1.15e -> %1.15e, %1.15e -> %1.15e, %1.15e -> %1.15e\n", xprev, x, yprev, y, zprev, z);
    //if(id == 0) printf("d:%1.15e, %1.15e, %1.15e\ns:%1.15e, %1.15e, %1.15e\n", x-xprev, y-yprev, z-zprev, xE, yE, zE);
    //if(id == 0) printf("%s %s %s\n", sign(x-xprev)==sign(xE)?"true ":"false",sign(y-yprev)==sign(yE)?"true ":"false",sign(z-zprev)==sign(zE)?"true ":"false");
    if(xprev == x && yprev == y && zprev == z) break;
  }
  x0[id] = xprev;
  y0[id] = yprev;
  z0[id] = zprev;
  x1[id] = x;
  y1[id] = y;
  z1[id] = z;
}

void kernel
tnp_k_implicit(global const float8 *a1, global const float8 *a2,                   // E, B coeff
      global float *x0, global float *y0, global float *z0, // initial pos
      global float *x1, global float *y1, global float *z1,
      const float Bcoeff, const float Ecoeff, // Bcoeff, Ecoeff
      const unsigned int n, const unsigned int ncalc // n, ncalc
      ) {
  uint id = get_global_id(0);
  uint prev_idx = UINT_MAX;
  float xprev = x0[id], yprev = y0[id], zprev = z0[id],
        x = x1[id], y = y1[id], z = z1[id];
  float8 temp, pos;
  float8 store0, store1, store2, store3, store4, store5;
  for (int t = 0; t < ncalc; t++) {
    if(x <= XLOW || x >= XHIGH || y <= YLOW || y >= YHIGH || z <= ZLOW || z >= ZHIGH) break;
    float xy = x * y, xz = x * z, yz = y * z, xyz = x * yz;
    uint idx = ((uint)((z - ZLOW) / DZ) * NZ + (uint)((y - YLOW) / DY)) * NY + (uint)((x - XLOW) / DX);
    // round down the cells - this is intentional
    idx *= 3;
    pos = (float8)(1.f, x, y, z, xy, xz, yz, xyz);

    // Is there no better way to do this? Why does float8 not have dot()?
    if (prev_idx != idx) {
      store0 = a1[idx]; store1 = a1[idx+1]; store2 = a1[idx+2]; store3 = a2[idx]; store4 = a2[idx+1]; store5 = a2[idx+2];
      prev_idx = idx;
    }
    temp = store0 * pos;
    float xE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store1 * pos;
    float yE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store2 * pos;
    float zE = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store3 * pos;
    float xP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store4 * pos;
    float yP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;
    temp = store5 * pos;
    float zP = temp.s0 + temp.s1 + temp.s2 + temp.s3 + temp.s4 + temp.s5 + temp.s6 + temp.s7;

    xP *= Bcoeff;
    yP *= Bcoeff;
    zP *= Bcoeff;
    xE *= Ecoeff;
    yE *= Ecoeff;
    zE *= Ecoeff;

    float xyP = xP * yP, yzP = yP * zP, xzP = xP * zP;
    float xxP = xP * xP, yyP = yP * yP, zzP = zP * zP;
    float b_det = 1.f / (1.f + xxP + yyP + zzP);

    float vx = (x - xprev); // / dt -> cancels out in the end
    float vy = (y - yprev);
    float vz = (z - zprev);

    xprev = x;
    yprev = y;
    zprev = z;

    //float vxnext = b_det * ((b_prime - 2 * (yyP + zzP)) * vx + (2 * (zP + xyP)) * vy + (2 * (xzP - yP)) * vz + (1.f + xxP) * xE + (zP + xyP) * yE + (xzP - yP) * zE);
    //float vynext = b_det * ((2 * (xyP - zP)) * vx + (b_prime - 2 * (xxP + zzP)) * vy + (2 * (xP + yzP)) * vz + (xyP - zP) * xE + (1.f + yyP) * yE + (yzP + xP) * zE);
    //float vznext = b_det * ((2 * (yP + xzP)) * vx + (2 * (yzP - xP)) * vy + (b_prime - 2 * (xxP + yyP)) * vz + (xzP + yP) * xE + (yzP - xP) * yE + (1.f + zzP) * zE);

    // Collect common terms
    // float vxxe = 2 * vx + xE,
    //       vyye = 2 * vy + yE,
    //       vzze = 2 * vz + zE;
    // float vxnext = vx + b_det * (-2 * vx * (yyP + zzP) + vyye * (zP + xyP) + vzze * (xzP - yP) + (1.f + xxP) * xE);
    // float vynext = vy + b_det * (vxxe * (xyP - zP) + -2 * vy * (xxP + zzP) + vzze * (xP + yzP) + (1.f + yyP) * yE);
    // float vznext = vz + b_det * (vxxe * (yP + xzP) + vyye * (yzP - xP) + -2 * vz * (xxP + yyP) + (1.f + zzP) * zE);

    // (v2 + v1) / 2 * dt, move the dt out
    // x = clamp(x + (vxnext + vx) / 2.f, XLOW, XHIGH);
    // y = clamp(y + (vynext + vy) / 2.f, YLOW, YHIGH);
    // z = clamp(z + (vznext + vz) / 2.f, ZLOW, ZHIGH);

    // Divide by 2 (move 0.5 for xE to Ecoeff)
    float vxxe = vx + xE,
          vyye = vy + yE,
          vzze = vz + zE;
    //x += vx + b_det * (-vx * (yyP + zzP) + vyye * (zP + xyP) + vzze * (xzP - yP) + (1.f + xxP) * xE);
    //y += vy + b_det * (vxxe * (xyP - zP) -  vy * (xxP + zzP) + vzze * (xP + yzP) + (1.f + yyP) * yE);
    //z += vz + b_det * (vxxe * (yP + xzP) + vyye * (yzP - xP) -  vz * (xxP + yyP) + (1.f * zzP) * zE);
    // do fma
    x += fma(b_det, fma(-vx, yyP + zzP, fma(vyye, zP + xyP, fma(vzze, xzP - yP, fma(xxP, xE, xE)))), vx);
    y += fma(b_det, fma(vxxe, xyP - xP, fma(-vy, xxP + zzP, fma(vzze, xP + yzP, fma(yyP, yE, yE)))), vy);
    z += fma(b_det, fma(vxxe, yP + xzP, fma(vyye, yzP - xP, fma(-vz, xxP + yyP, fma(zzP, zE, zE)))), vz);
  }
  x0[id] = xprev;
  y0[id] = yprev;
  z0[id] = zprev;
  x1[id] = x;
  y1[id] = y;
  z1[id] = z;
}