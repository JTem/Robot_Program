
class InverseKinematik {

  Jama.Matrix[] T0 = new Jama.Matrix[7];
  Jama.Matrix Jb = Jama.Matrix.identity(4,4);

  //float q1=radians(1.1), q2=radians(1.2), q3=radians(1.3), q4=radians(1.4), q5=radians(1.5), q6=radians(1.6);
  float q1=0, q2=0, q3=0, q4=0, q5=0.1, q6=0;
  Jama.Matrix t11, t12, t13, t14, t15, t16, t66, t65, t64, t63, t62;
  float theta_MatrixLog;




  public InverseKinematik() {
    body_forwardKinematik(q1, q2, q3, q4, q5, q6);
  }

  //antisymmetric / skew symmetric Matrix
  public Jama.Matrix asMatrix(double x1, double x2, double x3) {

    double[][] asm={
      {0, -x3, x2}, 
      {x3, 0, -x1}, 
      {-x2, x1, 0}};

    Jama.Matrix ASM = new Jama.Matrix(asm);
    return ASM;
  }



  public Jama.Matrix asMatrix(Jama.Matrix w) {

    double[][] asm={
      {0, -w.get(2, 0), w.get(1, 0)}, 
      {w.get(2, 0), 0, -w.get(0, 0)}, 
      {-w.get(1, 0), w.get(0, 0), 0}};

    Jama.Matrix ASM = new Jama.Matrix(asm);
    return ASM;
  }


  public Jama.Matrix rodriguesFormula(double x1, double x2, double x3, float theta) {

    Jama.Matrix I = Jama.Matrix.identity(3, 3);
    Jama.Matrix asm = asMatrix(x1, x2, x3);
    Jama.Matrix asm_square = asm.times(asm);
    Jama.Matrix m1 = asm.times(sin(theta));
    Jama.Matrix m2 = asm_square.times((1-cos(theta)));
    Jama.Matrix rod = I.plus(m1.plus(m2));
    return rod;
  }




  public Jama.Matrix matrixExpo(Jama.Matrix s, float theta) {

    Jama.Matrix I = Jama.Matrix.identity(3, 3);
    Jama.Matrix R = rodriguesFormula(s.get(0, 0), s.get(1, 0), s.get(2, 0), theta);
    Jama.Matrix I_theta = I.times(theta);
    Jama.Matrix asm = asMatrix(s.get(0, 0), s.get(1, 0), s.get(2, 0));
    Jama.Matrix asm_square = asm.times(asm);
    Jama.Matrix m1 = asm.times((1-cos(theta)));
    Jama.Matrix m2 = asm_square.times((theta-sin(theta)));
    Jama.Matrix m3 = I_theta.plus(m1);
    Jama.Matrix m4 = m3.plus(m2);
    double[][] sv = {
      {s.get(3, 0)}, 
      {s.get(4, 0)}, 
      {s.get(5, 0)}};

    Jama.Matrix Sv = new Jama.Matrix(sv); 
    Jama.Matrix P = m4.times(Sv);

    double [][] expo = {
      {R.get(0, 0), R.get(0, 1), R.get(0, 2), P.get(0, 0)}, 
      {R.get(1, 0), R.get(1, 1), R.get(1, 2), P.get(1, 0)}, 
      {R.get(2, 0), R.get(2, 1), R.get(2, 2), P.get(2, 0)}, 
      {0, 0, 0, 1}};
    Jama.Matrix Expo = new Jama.Matrix(expo);
    return Expo;
  }




  public Jama.Matrix matrixLog(Jama.Matrix T) {

    Jama.Matrix R = T.getMatrix(0, 2, 0, 2);
    Jama.Matrix RT = R.transpose();
    Jama.Matrix RRT = R.minus(RT);
    Jama.Matrix P = T.getMatrix(0, 2, 3, 3);
    double p_norm = P.norm2();
    Jama.Matrix v;
    Jama.Matrix I = Jama.Matrix.identity(3, 3);
    Jama.Matrix RI = R.minus(I);
    double trRi = RI.trace();
    Jama.Matrix asmW = Jama.Matrix.identity(3, 3);
    double trR = R.trace();
    double theta;

    if (abs((float)trRi) <= 0.001) {

      double[][] w0 = {{0}, {0}, {0}};
      Jama.Matrix W0 = new Jama.Matrix(w0);
      asmW = asMatrix(W0);

      double factorV = 1/p_norm;
      v = P.times(factorV);
      theta = p_norm;
    } else if (trR + 1 <= 0.001) {
      theta = PI;

      if (R.get(2, 2)>=-0.98) {

        float factor1 = 1 / sqrt(2*(1+(float)R.get(2, 2)));
        double[][] w1 = {{R.get(0, 2)}, {R.get(1, 2)}, {1+R.get(2, 2)}};
        Jama.Matrix W1 = new Jama.Matrix(w1);
        W1.timesEquals(factor1);
        asmW = asMatrix(W1);
      }
      if (R.get(1, 1)>=-0.98) {

        float factor2 = 1 / sqrt(2*(1+(float)R.get(1, 1)));
        double[][] w2 = {{R.get(0, 1)}, {1+R.get(1, 1)}, {R.get(2, 1)}};
        Jama.Matrix W2 = new Jama.Matrix(w2);
        W2.timesEquals(factor2);
        asmW = asMatrix(W2);
      }
      if (R.get(0, 0)>=-0.98) {

        float factor3 = 1 / sqrt(2*(1+(float)R.get(0, 0)));
        double[][] w3 = {{1+R.get(0, 0)}, {R.get(1, 0)}, {R.get(2, 0)}};
        Jama.Matrix W3 = new Jama.Matrix(w3);
        W3.timesEquals(factor3);
        asmW = asMatrix(W3);
      }

      double factorT = 1/theta;
      double factorasm = -0.5;
      float cot = 1/tan(0.5*(float)theta);
      double factorasm_sq = factorT-0.5*cot;
      Jama.Matrix asm_sq = asmW.times(asmW);
      Jama.Matrix I_T = I.times(factorT);
      Jama.Matrix m1 = asmW.times(factorasm);
      Jama.Matrix m2 = asm_sq.times(factorasm_sq);
      Jama.Matrix m3 = I_T.plus(m1);
      Jama.Matrix m4 = m3.plus(m2);

      v=m4.times(P);
    } else {

      theta = acos(0.5*((float)trR-1));
      float factor = 1 / (2*sin((float)theta));
      asmW = RRT.times(factor);

      double factorT = 1/theta;
      double factorasm = -0.5;
      float cot = 1/tan(0.5*(float)theta);
      double factorasm_sq = factorT-0.5*cot;
      Jama.Matrix asm_sq = asmW.times(asmW);
      Jama.Matrix I_T = I.times(factorT);
      Jama.Matrix m1 = asmW.times(factorasm);
      Jama.Matrix m2 = asm_sq.times(factorasm_sq);
      Jama.Matrix m3 = I_T.plus(m1);
      Jama.Matrix m4 = m3.plus(m2);

      v=m4.times(P);
    }


    double[][] twist = {
      {asmW.get(0, 0), asmW.get(0, 1), asmW.get(0, 2), v.get(0, 0)}, 
      {asmW.get(1, 0), asmW.get(1, 1), asmW.get(1, 2), v.get(1, 0)}, 
      {asmW.get(2, 0), asmW.get(2, 1), asmW.get(2, 2), v.get(2, 0)}, 
      {0, 0, 0, 0}};

    Jama.Matrix V = new Jama.Matrix(twist);
    V.timesEquals(theta);
    theta_MatrixLog = (float)theta;
    return V;
  }


  public Jama.Matrix adJoint(Jama.Matrix T) {

    Jama.Matrix asm = asMatrix(T.get(0, 3), T.get(1, 3), T.get(2, 3));
    Jama.Matrix R = T.getMatrix(0, 2, 0, 2);
    Jama.Matrix pR = asm.times(R);

    double[][] adj = {
      {R.get(0, 0), R.get(0, 1), R.get(0, 2), 0, 0, 0}, 
      {R.get(1, 0), R.get(1, 1), R.get(1, 2), 0, 0, 0}, 
      {R.get(2, 0), R.get(2, 1), R.get(2, 2), 0, 0, 0}, 
      {pR.get(0, 0), pR.get(0, 1), pR.get(0, 2), R.get(0, 0), R.get(0, 1), R.get(0, 2)}, 
      {pR.get(1, 0), pR.get(1, 1), pR.get(1, 2), R.get(1, 0), R.get(1, 1), R.get(1, 2)}, 
      {pR.get(2, 0), pR.get(2, 1), pR.get(2, 2), R.get(2, 0), R.get(2, 1), R.get(2, 2)}};

    Jama.Matrix Adj = new Jama.Matrix(adj);

    return Adj;
  }


  void body_forwardKinematik(double _q1, double _q2, double _q3, double _q4, double _q5, double _q6) {

    t11 = matrixExpo(Parameter.B1, (float)_q1);
    t12 = t11.times(matrixExpo(Parameter.B2, (float)_q2));
    t13 = t12.times(matrixExpo(Parameter.B3, (float)_q3));
    t14 = t13.times(matrixExpo(Parameter.B4, (float)_q4));
    t15 = t14.times(matrixExpo(Parameter.B5, (float)_q5));
    t16 = t15.times(matrixExpo(Parameter.B6, (float)_q6));

    t66 = matrixExpo(Parameter.B6, -(float)_q6);
    t65 = t66.times(matrixExpo(Parameter.B5, -(float)_q5));
    t64 = t65.times(matrixExpo(Parameter.B4, -(float)_q4));
    t63 = t64.times(matrixExpo(Parameter.B3, -(float)_q3));
    t62 = t63.times(matrixExpo(Parameter.B2, -(float)_q2));

    T0[1] = t11;
    T0[2] = Parameter.M2.times(t12);
    T0[3] = Parameter.M3.times(t13);
    T0[4] = Parameter.M4.times(t14);
    T0[5] = Parameter.M5.times(t15);
    T0[6] = Parameter.M6.times(t16);
  }


  public Jama.Matrix get_current_Trans() {
    return T0[6];
  }

  void newton_raphson_ik(Jama.Matrix T, Jama.Matrix _qguess) {

    Jama.Matrix Tbd = T_inv(T0[6]).times(T);
    Jama.Matrix Twist = matrixLog(Tbd);
    if (_qguess.get(3, 0) > PI || _qguess.get(3, 0)<-PI)
      _qguess.set(3, 0, 0);
    Jama.Matrix q = _qguess;
    double eps_w = 0.01;  
    double eps_v = 0.5;
    int i = 0;
    double[][] t = {
      {Twist.get(2, 1)}, 
      {Twist.get(0, 2)}, 
      {Twist.get(1, 0)}, 
      {Twist.get(0, 3)}, 
      {Twist.get(1, 3)}, 
      {Twist.get(2, 3)}};

    Jama.Matrix Vb = new Jama.Matrix(t);
    Jama.Matrix wb = Vb.getMatrix(0, 2, 0, 0);
    Jama.Matrix vb = Vb.getMatrix(3, 5, 0, 0);

    while (vb.norm2() > eps_v || wb.norm2() > eps_w) {

      Jama.Matrix dq = pseudoinverse(body_Jacobian()).times(Vb);

      q.plusEquals(dq);
      q = set_Q_boundaries(q);
      set_q(q);
      body_forwardKinematik(q1, q2, q3, q4, q5, q6);

      Tbd = T_inv(T0[6]).times(T);
      Twist = matrixLog(Tbd);

      double[][] twist = {
        {Twist.get(2, 1)}, 
        {Twist.get(0, 2)}, 
        {Twist.get(1, 0)}, 
        {Twist.get(0, 3)}, 
        {Twist.get(1, 3)}, 
        {Twist.get(2, 3)}};

      Vb = new Jama.Matrix(twist);
      wb = Vb.getMatrix(0, 2, 0, 0);
      vb = Vb.getMatrix(3, 5, 0, 0);

      i++;
      if (i>20) break;
    }

    i=0;
    q = set_Q_boundaries(q);

    set_q(q);
  }


  void set_q(Jama.Matrix _q) {

    q1 = (float)_q.get(0, 0);
    q2 = (float)_q.get(1, 0);
    q3 = (float)_q.get(2, 0);
    q4 = (float)_q.get(3, 0);
    q5 = (float)_q.get(4, 0);
    q6 = (float)_q.get(5, 0);
  }

  public Jama.Matrix set_Q_boundaries(Jama.Matrix _q) {

    double _q1 = _q.get(0, 0);
    double _q2 = _q.get(1, 0);
    double _q3 = _q.get(2, 0);
    double _q4 = _q.get(3, 0);
    double _q5 = _q.get(4, 0);
    double _q6 = _q.get(5, 0);

    if (_q1 > Parameter.q1_pos) _q1=Parameter.q1_pos;
    if (_q1 < Parameter.q1_neg) _q1=Parameter.q1_neg;

    if (_q2 > Parameter.q2_pos) _q2=Parameter.q2_pos;
    if (_q2 < Parameter.q2_neg) _q2=Parameter.q2_neg;

    if (_q3 > Parameter.q3_pos) _q3=Parameter.q3_pos;
    if (_q3 < Parameter.q3_neg) _q3=Parameter.q3_neg;

    if (_q4 > Parameter.q4_pos) _q4=Parameter.q4_pos;
    if (_q4 < Parameter.q4_neg) _q4=Parameter.q4_neg;

    if (_q5 > Parameter.q5_pos) _q5=Parameter.q5_pos;
    if (_q5 < Parameter.q5_neg) _q5=Parameter.q5_neg;

    double[][] r = {{_q1}, {_q2}, {_q3}, {_q4}, {_q5}, {_q6}};
    Jama.Matrix res = new Jama.Matrix(r);
    return res;
  }


  public Jama.Matrix get_qVec() {
    double[][] q = {
      {q1}, 
      {q2}, 
      {q3}, 
      {q4}, 
      {q5}, 
      {q6}};
    Jama.Matrix Q = new Jama.Matrix(q);

    return Q;
  }


  public Jama.Matrix body_Jacobian() {

    Jama.Matrix Jb1 = adJoint(t62).times(Parameter.B1);
    Jama.Matrix Jb2 = adJoint(t63).times(Parameter.B2);
    Jama.Matrix Jb3 = adJoint(t64).times(Parameter.B3);
    Jama.Matrix Jb4 = adJoint(t65).times(Parameter.B4);
    Jama.Matrix Jb5 = adJoint(t66).times(Parameter.B5);
    Jama.Matrix Jb6 = Parameter.B6;


    double[][] jb ={
      {Jb1.get(0, 0), Jb2.get(0, 0), Jb3.get(0, 0), Jb4.get(0, 0), Jb5.get(0, 0), Jb6.get(0, 0)}, 
      {Jb1.get(1, 0), Jb2.get(1, 0), Jb3.get(1, 0), Jb4.get(1, 0), Jb5.get(1, 0), Jb6.get(1, 0)}, 
      {Jb1.get(2, 0), Jb2.get(2, 0), Jb3.get(2, 0), Jb4.get(2, 0), Jb5.get(2, 0), Jb6.get(2, 0)}, 
      {Jb1.get(3, 0), Jb2.get(3, 0), Jb3.get(3, 0), Jb4.get(3, 0), Jb5.get(3, 0), Jb6.get(3, 0)}, 
      {Jb1.get(4, 0), Jb2.get(4, 0), Jb3.get(4, 0), Jb4.get(4, 0), Jb5.get(4, 0), Jb6.get(4, 0)}, 
      {Jb1.get(5, 0), Jb2.get(5, 0), Jb3.get(5, 0), Jb4.get(5, 0), Jb5.get(5, 0), Jb6.get(5, 0)}};

    Jb = new Jama.Matrix(jb);

    return Jb;
  }





  public Jama.Matrix pseudoinverse(Jama.Matrix _J) {
    Jama.Matrix I = Jama.Matrix.identity(_J.getRowDimension(), _J.getColumnDimension());
    double xi = 0.1;
    Jama.Matrix JT = _J.transpose();
    Jama.Matrix JJT = JT.times(_J);
    JJT.plusEquals(I.times(xi));
    Jama.Matrix JJT_inv = JJT.inverse();
    Jama.Matrix p = JJT_inv.times(JT);
    return p;
  }



  public Jama.Matrix T_inv(Jama.Matrix T) {

    Jama.Matrix R = T.getMatrix(0, 2, 0, 2);
    Jama.Matrix RT= R.transpose();
    Jama.Matrix P = T.getMatrix(0, 2, 3, 3);
    Jama.Matrix RTP = RT.times(P);

    double[][] ti = {
      {RT.get(0, 0), RT.get(0, 1), RT.get(0, 2), -RTP.get(0, 0)}, 
      {RT.get(1, 0), RT.get(1, 1), RT.get(1, 2), -RTP.get(1, 0)}, 
      {RT.get(2, 0), RT.get(2, 1), RT.get(2, 2), -RTP.get(2, 0)}, 
      {0, 0, 0, 1}};

    Jama.Matrix Ti = new Jama.Matrix(ti);

    return Ti;
  }



  void printRMatrix(Jama.Matrix m) {

    println ((float)m.get(0, 0)+"  "+(float)m.get(0, 1)+"  "+(float)m.get(0, 2)+"  "+(float)m.get(0, 3));
    println ((float)m.get(1, 0)+"  "+(float)m.get(1, 1)+"  "+(float)m.get(1, 2)+"  "+(float)m.get(1, 3));
    println ((float)m.get(2, 0)+"  "+(float)m.get(2, 1)+"  "+(float)m.get(2, 2)+"  "+(float)m.get(2, 3));
    println("------------------------------------------------------------------------");
  } 



  void showKS(int jointNumber) {

    double[][] ks = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    Jama.Matrix getKS = new Jama.Matrix(ks);

    switch(jointNumber) {

    case 0:
      break;

    case 1:
      getKS = T0[1];
      break;

    case 2:
      getKS = T0[2];
      break;

    case 3:
      getKS = T0[3];
      break;

    case 4:
      getKS = T0[4];
      break;

    case 5:
      getKS = T0[5];
      break;

    case 6:
      getKS = T0[6];
      break;
    }

    float mx=(float)getKS.get(0, 0);
    float my=(float)getKS.get(1, 0);
    float mz=(float)getKS.get(2, 0);
    float lx=(float)getKS.get(0, 1);
    float ly=(float)getKS.get(1, 1);
    float lz=(float)getKS.get(2, 1);
    float nx=(float)getKS.get(0, 2);
    float ny=(float)getKS.get(1, 2);
    float nz=(float)getKS.get(2, 2);
    float x =(float)getKS.get(0, 3);
    float y =(float)getKS.get(1, 3);
    float z =(float)getKS.get(2, 3);

    float size=80;

    pushMatrix();
    strokeWeight(1);
    stroke(255, 0, 0);
    line(x, y, z, x+mx*size, y+my*size, z+ mz*size);   
    stroke(0, 255, 0); 
    line(x, y, z, x+lx*size, y+ly*size, z+lz*size);
    stroke(0, 0, 255);
    line(x, y, z, x+nx*size, y+ny*size, z+nz*size);
    popMatrix();
  }
}
