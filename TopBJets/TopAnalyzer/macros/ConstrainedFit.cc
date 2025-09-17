#include <iostream>
#include <assert.h>

using namespace std;

#include "TVectorD.h"
#include "TArrayD.h"
#include "ConstrainedFit.hh"

#define FIT_VERBOSE 0

ConstrainedFit::ConstrainedFit(int n,int m,int r) {
  _ncon = r;
  _npar = m;
  _maxtrk = n;
  _ntrk = 0;
  _order = 0;
  _status = 1;
  _verbose = FIT_VERBOSE;
  _abz = ABZ;
  _par = TVectorD();//'Error in <TVectorT<double>::Allocate>: nrows=-5',_par = TVectorD(m,0) to _par = TVectorD();
  _par.ResizeTo(m);
  _cov = TMatrixD();
  _cov.ResizeTo(m,m);
}

ConstrainedFit::~ConstrainedFit() { }

void ConstrainedFit::AddData(const TVectorD &v,const TMatrixD &c,
                             double m,unsigned int f) {
  if(_ntrk < _maxtrk) {
    _mask.push_back(f);
    _mass.push_back(m);
    _xx.push_back(v);
    _vv.push_back(c);
    _x.push_back(v);
    _ntrk += 1;
    _order += v.GetNrows();
  }
  else {
    cout<<"ConstrainedFit::AddData() - Too many tracks, maxtrk="
         <<_maxtrk<<"."<<endl;
    _status = -1;
  }
}

void ConstrainedFit::AddTrack(const double *p,const double *c,
                              double m,unsigned int f) {
  assert(p!=NULL&&c!=NULL);
  TVectorD par(5);
  TMatrixD cov(5,5);
  int k = 0;
  for(int i=0; i<5; i++) {
    par[i] = p[i];
    for(int j=0; j<5; j++) {
      cov[i][j] = c[k++];
    }
  }
  AddData(par,cov,m,f);
}

void ConstrainedFit::GetVertexMass(unsigned int iv,double &m,double &em) const {
//
//  Vertices made from tracks + pseudo-tracks need their own code to
//  calculate this.  Otherwise the base class is sufficient.
//
  double p[4] = { 0, 0, 0, 0 };
  for ( int i=0; i<NTrack(); i++ ) {
    if ( Vertex(i) == iv ) {
      double pt = fabs(ABz()/_x[i][0]);
      double px = pt*cos(_x[i][1]);
      double py = pt*sin(_x[i][1]);
      double pz = pt*_x[i][3];
      double m = Mass(i);
      double e = sqrt(pt*pt+pz*pz+m*m);
      p[0] += px;
      p[1] += py;
      p[2] += pz;
      p[3] += e;
    }
  }
  m = sqrt(p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2]);
  // We won't calculate the uncertainty yet...
  em = 0;
}

double ConstrainedFit::VertexMass(unsigned int iv) const {
  double m, em;
  GetVertexMass(iv,m,em);
  return m;
}

TVectorD ConstrainedFit::dGdP(int ipar) {
  const double h[3] = {0.001, 0.0005, 0.0001};
  TMatrixDSym hh(3);
  hh[0][0] = 3;
  hh[1][0] = h[0]+h[1]+h[2];
  hh[1][1] = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  hh[2][0] = hh[1][1];
  hh[2][1] = pow(h[0],3)+pow(h[1],3)+pow(h[2],3);
  hh[2][2] = pow(h[0],4)+pow(h[1],4)+pow(h[2],4);
  double ierr;
  hh.Invert(&ierr);
  if(ierr == 0) {
    cerr << "ConstrainedFit::dGdP() - Error inverting hh." << endl;
    assert(0);
  }
  TVectorD dg[3];
  TVectorD newpar = _par;
  for(int i=0; i<3; i++) {
    newpar[ipar] = _par[ipar] - h[i]/2;
    TVectorD g0 = G(_x,newpar);
    newpar[ipar] = _par[ipar] + h[i]/2;
    TVectorD g1 = G(_x,newpar);
    dg[i].ResizeTo(NConst()); // Added for TVectorD
    dg[i] = (g1-g0)*(1/h[i]);
  }
  TVectorD a(NConst()); // NConst() = 7
  for ( int i=0; i<NConst(); i++ ) {
    double yh[3];
    yh[0] = dg[0][i] + dg[1][i] + dg[2][i];
    yh[1] = dg[0][i]*h[0] + dg[1][i]*h[1] + dg[2][i]*h[2];
    yh[2] = dg[0][i]*h[0]*h[0] + dg[1][i]*h[1]*h[1] + dg[2][i]*h[2]*h[2];
    a[i] = hh[0][0]*yh[0] + hh[0][1]*yh[1] + hh[0][2]*yh[2];
  }
  return a;
}

TVectorD ConstrainedFit::dGdX(int itrk,int ipar) {
  const double h[3] = {0.0001, 0.00005, 0.00001};
  TMatrixDSym hh(3);
  hh[0][0] = 3;
  hh[1][0] = h[0]+h[1]+h[2];
  hh[1][1] = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  hh[2][0] = hh[1][1];
  hh[2][1] = pow(h[0],3)+pow(h[1],3)+pow(h[2],3);
  hh[2][2] = pow(h[0],4)+pow(h[1],4)+pow(h[2],4);
  double ierr;
  hh.Invert(&ierr);
  if(ierr == 0) {
    cerr<<"ConstrainedFit::dGdX() - Error inverting hh."<<endl;
    assert(0);
  }
  TVectorD dg[3];
  std::vector<TVectorD> newx = _x;
  for ( int i=0; i<3; i++ ) {
    newx[itrk][ipar] = _x[itrk][ipar] - h[i]/2;
    TVectorD g0 = G(newx,_par);
    newx[itrk][ipar] = _x[itrk][ipar] + h[i]/2;
    TVectorD g1 = G(newx,_par);
    dg[i].ResizeTo(NConst()); // Added for TVectorD
    dg[i] = (g1-g0)*(1/h[i]);
  }
  TVectorD a(NConst());
  for ( int i=0; i<NConst(); i++ ) {
    double yh[3];
    yh[0] = dg[0][i] + dg[1][i] + dg[2][i];
    yh[1] = dg[0][i]*h[0] + dg[1][i]*h[1] + dg[2][i]*h[2];
    yh[2] = dg[0][i]*h[0]*h[0] + dg[1][i]*h[1]*h[1] + dg[2][i]*h[2]*h[2];
    a[i] = hh[0][0]*yh[0] + hh[0][1]*yh[1] + hh[0][2]*yh[2];
  }
  return a;
}

int ConstrainedFit::CheckB() {
  int bwarn = 0;
  const double eps = 0.001;
  TMatrixD b = B(_x,_par);
  for(int i=0; i<NPar(); i++) { // NPar() = 6
    TVectorD dgdp = dGdP(i);
    for(int j=0; j<NConst(); j++) {
      if(_verbose >= 2 || (fabs(dgdp[j]-b[j][i]) > eps*fabs(dgdp[j]) &&
                            fabs(dgdp[j]-b[j][i]) > eps*fabs(b[j][i]))) {
        //cout<<"ConstrainedFit::B["<<j<<","<<i<<"] = "
             //<<b[j][i]<<"  dGdP = "<<dgdp[j]<<endl;
        if(fabs(dgdp[j]-b[j][i]) > eps) bwarn += 1;
      }
    }
  }
  return bwarn;
}

int ConstrainedFit::CheckA() {
  int awarn = 0;
  const double eps = 0.001;
  TMatrixD a = A(_x,_par);
  int l = 0;
  for(int i=0; i<NTrack(); i++) {
    for( int j=0; j<_x[i].GetNrows(); j++) {
      TVectorD dgdx = dGdX(i,j);
      for(int k=0; k<NConst(); k++) {
        if(_verbose >= 2 || fabs(dgdx[k]-a[l+j][k]) > eps*fabs(dgdx[k])) {
          //cout<<"ConstrainedFit::A["<<l+j<<"]["<<k
              //<<"] = "<<a[l+j][k]<<"  dG("<<k<<")/dx("
              //<<i<<","<<j<<") = "<<dgdx[k]<<endl;
          if(fabs(dgdx[k]-a[l+j][k]) > eps*fabs(dgdx[k])) awarn += 1;
        }
      }
    }
    l += _x[i].GetNrows();
  }
  return awarn;
}

int ConstrainedFit::InitialEstimates(const std::vector<TVectorD> &x, TVectorD &par) {
  for(int i=0; i<NPar(); i++) {
    par[i] = 0;
  }
  return 1;
}

void ConstrainedFit::CalculateDerivatives(const std::vector<TVectorD> &x,
                                                                  const TVectorD &par) { }

int ConstrainedFit::Fit(int niter) {
  int ierr = InitialEstimates(_x,_par);
  if(ierr != 0) {
    cout <<
    "ConstrainedFit::Fit() - Error calculating initial estimates."<<endl;
    return ierr;
  }

  if(_verbose >= 2) {
    cout <<"ConstrainedFit::Fit() - "<<_ntrk<<" tracks, "
                                     <<_ncon<<" constraints, " 
                                     <<_npar<<" parameters."<<endl;
    cout<<"Order = " << _order << endl;
    cout<<"niter = "<<niter<<endl;
    //_x.front().Print(); _par.Print();
    TVectorD p1 = _xx.at(0); TVectorD p2 = _xx.at(1); p1.Print(); p2.Print(); 
    if(_xx.size() == 3) {TVectorD p3 = _xx.at(2); p3.Print();}
    TMatrixD c1 = _vv.at(0); TMatrixD c2 = _vv.at(1); c1.Print(); c2.Print();
    if(_xx.size() == 3) {TMatrixD c3 = _vv.at(2); c3.Print();}
  }

  if(_verbose && (ierr = CheckB()) > 0 ) {
    cout<<"ConstrainedFit::Fit() - Derivative matrix B: "<<ierr
         <<" warnings."<<endl;
  }
  else if(_verbose == 2) {
    cout<<"ConstrainedFit::Fit() - Derivative matrix B: no warnings."<<endl;
  }

  if(_verbose && (ierr = CheckA()) > 0) {
    cout<<"ConstrainedFit::Fit() - Derivative matrix A: "<<ierr
         <<" warnings." <<endl;
  }
  else if(_verbose == 2) {
    cout<<"ConstrainedFit::Fit() - Derivative matrix A: no warnings."<<endl;
  }

  TMatrixD v(_order,_order);
  int l = 0;
  for(int i=0; i<NTrack(); i++) {
    _x[i] = _xx[i];
    for(int j=0; j<_xx[i].GetNrows(); j++) {
      for(int k=0; k<_xx[i].GetNrows(); k++) {
        v[l+j][l+k] = _vv[i][j][k];
      }
    }
    l += _xx[i].GetNrows();
  }

  int iteration = 0;
  TMatrixD w(NConst(),NConst()), bwb(NPar(),NPar());
  TMatrixD a, b;
  TVectorD d(NPar()), e(5*NTrack());
  do {
    CalculateDerivatives(_x,_par);
    a.ResizeTo(5*NTrack(),NConst());
    a = A(_x,_par);
    TMatrixD at(TMatrixD::kTransposed,a);
    w = at*v*a;
    //w.Print(); at.Print(); v.Print(); a.Print();
    double ierr;
    w.Invert(&ierr); // This changes w to w^-1
    if(ierr == 0) {
      cout<<"ConstrainedFit::Fit() - Error inverting matrix W."<<endl;
      _status = -1;
      return _status;
    }
    b.ResizeTo(NConst(),NPar());
    b = B(_x,_par);
    TMatrixD bt(TMatrixD::kTransposed,b);
    bwb = bt*w*b;
    //bwb.Print(); bt.Print(); w.Print(); b.Print();
    bwb.Invert(&ierr);
    if(ierr == 0) {
      cout<<"ConstrainedFit::Fit() - Error inverting matrix BWB."<<endl;
      _status = -1;
      return _status;
    }
    TVectorD g = G(_x,_par);
    if(iteration > 0) {
      TVectorD dx(_order);
      int l = 0;
      for(int i=0; i<NTrack(); i++) {
        for(int j=0; j<_x[i].GetNrows(); j++) {
          dx[l+j] = _xx[i][j] - _x[i][j];
        }
        l += _x[i].GetNrows();
      }
      g += at*dx;
      //at.Print(); dx.Print();
    }
    //g.Print();
    d = -1.0*bwb*bt*w*g; // (6x6 * 6x7 * 7x7 * 7x1) -> (6x1)
    e = v*a*w*(g+b*d); // (10x1)
    //d.Print(); e.Print();
    for(int i=0; i<NPar(); i++) {
      _par[i] += d[i]; // Check
    }
    int l = 0;
    for(int i=0; i<NTrack(); i++) {
      for(int j=0; j<_x[i].GetNrows(); j++) {
        _x[i][j] -= e[l+j]; // Check
      }
      l += _x[i].GetNrows();
    }
  }
  while (++iteration < niter);
  _cov = bwb;
  _chi2 = 0;
  l = 0;
  for(int i=0; i<NTrack(); i++) {
    TMatrixD v = _vv[i];
    double ierr;
    v.Invert(&ierr);
    if(ierr == 0) {
      cout<<"ConstrainedFit::Fit() - Error inverting matrix V("
             <<i<<")."<<endl;
      _status = -1;
      return _status;
    }
    TVectorD y = e.GetSub(l,l+_x[i].GetNrows()-1); // Check
    TArrayD data(5);
    for(int i=0; i<y.GetNrows(); i++) {
      data[i] = y[i];
    }
    TMatrixD dx; dx.Use(5,1,data.GetArray()); TMatrixD dxt(TMatrixD::kTransposed,dx);
    TMatrixD x = dxt*v*dx;
    //x.Print(); dxt.Print(); v.Print(); dx.Print();
    _chi2 += x[0][0];
    l += _x[i].GetNrows();
  }
  _status = 0;
  return _status;
}

Fit3D::Fit3D(int n) : ConstrainedFit(n,6,3+2*n), mass_constraint(0) { }

int Fit3D::InitialEstimates(const std::vector<TVectorD> &x,TVectorD &par) {
  par[0] = 0.1;
  par[1] = 0;
  par[2] = 0;
  par[3] = 0;
  par[4] = 0;
  par[5] = 0;
  for(int i=0; i<NTrack(); i++) {
    double pt = fabs(ABz()/x[i][0]);
    double px = pt*cos(x[i][1]);
    double py = pt*sin(x[i][1]);
    double pz = pt*x[i][3];
    par[2] += x[i][4];
    par[3] += px;
    par[4] += py;
    par[5] += pz;
  }
  par[2] /= NTrack();
  return 0;
}

TVectorD Fit3D::G(const std::vector<TVectorD> &x,
                              const TVectorD &p) const {
  TVectorD g(NConst());
  double s = p[0];
  double d0 = p[1];
  double z0 = p[2];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  double xv = (s*qx-d0*qy)/qxy;
  double yv = (s*qy+d0*qx)/qxy;
  double zv = z0 + qz*s/qxy;
  g[0] = -qx;
  g[1] = -qy;
  g[2] = -qz;
  for ( int i=0; i<NTrack(); i++ ) {
    double pt = fabs(ABz()/x[i][0]);
    double px = pt*cos(x[i][1]);
    double py = pt*sin(x[i][1]);
    double pz = pt*x[i][3];
    g[0] += px;
    g[1] += py;
    g[2] += pz;
    double s0 = xv*cos(x[i][1]) + yv*sin(x[i][1]);
    g[3+2*i] = x[i][2] + xv*sin(x[i][1]) - yv*cos(x[i][1]);
    g[3+2*i+1] = x[i][4] - zv + s0*x[i][3];
  }
  return g;
}

TMatrixD Fit3D::A(const std::vector<TVectorD> &x,
                             const TVectorD &p) const {
  double s = p[0];
  double d0 = p[1];
  double z0 = p[2];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  double xv = (s*qx-d0*qy)/qxy;
  double yv = (s*qy+d0*qx)/qxy;
  double zv = z0 + s*qz/qxy;
  TMatrixD a(5*NTrack(),NConst());
  for(int i=0; i<NTrack(); i++) {
    double pt = fabs(ABz()/x[i][0]);
    double px = pt*cos(x[i][1]);
    double py = pt*sin(x[i][1]);
    double pz = pt*x[i][3];
    double s0 = xv*cos(x[i][1]) + yv*sin(x[i][1]);
    a[5*i][0] = -px/x[i][0];
    a[5*i+1][0] = -py;

    a[5*i][1] = -py/x[i][0];
    a[5*i+1][1] = px;

    a[5*i][2] = -pz/x[i][0];
    a[5*i+3][2] = pz/x[i][3];

    a[5*i+1][3+2*i] = xv*cos(x[i][1]) + yv*sin(x[i][1]);
    a[5*i+2][3+2*i] = 1.0;

    a[5*i+1][3+2*i+1] = (-xv*sin(x[i][1])+yv*cos(x[i][1]))*x[i][3];
    a[5*i+3][3+2*i+1] = s0;
    a[5*i+4][3+2*i+1] = 1.0;
  }
  return a;
}

TMatrixD Fit3D::B(const std::vector<TVectorD> &x,
                             const TVectorD &p) const {
  double s = p[0];
  double d0 = p[1];
  double z0 = p[2];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  TMatrixD b(NConst(),6);
  b[0][3] = -1.0;
  b[1][4] = -1.0;
  b[2][5] = -1.0;
  for(int i=0; i<NTrack(); i++) {
    b[3+2*i][0] = (qx*sin(x[i][1])-qy*cos(x[i][1]))/qxy;
    b[3+2*i][1] = (-qy*sin(x[i][1])-qx*cos(x[i][1]))/qxy;
    b[3+2*i][3] = qy*((s*qy+d0*qx)*sin(x[i][1])/qxy+
                      (s*qx-d0*qy)*cos(x[i][1])/qxy)/(qxy*qxy);
    b[3+2*i][4] = -qx*((s*qy+d0*qx)*sin(x[i][1])/qxy+
                       (s*qx-d0*qy)*cos(x[i][1])/qxy)/(qxy*qxy);
    b[3+2*i+1][0] = x[i][3]*(qx*cos(x[i][1])+qy*sin(x[i][1]))/qxy - qz/qxy;
    b[3+2*i+1][1] = x[i][3]*(qx*sin(x[i][1])-qy*cos(x[i][1]))/qxy;
    b[3+2*i+1][2] = -1;
    b[3+2*i+1][3] = x[i][3]*qy*((s*qy+d0*qx)*cos(x[i][1])/qxy-
                                (s*qx-d0*qy)*sin(x[i][1])/qxy)/(qxy*qxy) +
                    s*qx*qz/pow(qxy,3);
    b[3+2*i+1][4] = -x[i][3]*qx*((s*qy+d0*qx)*cos(x[i][1])/qxy-
                                 (s*qx-d0*qy)*sin(x[i][1])/qxy)/(qxy*qxy) +
                    s*qy*qz/pow(qxy,3);
    b[3+2*i+1][5] = -s/qxy;
  }
  return b;
}

// ---------------------------------------------------------------------
//    Single vertex fit in 3D for CMS track parameterization
//    M Jones - 6-Apr-2022
//  ---------------------------------------------------------------------
CmsFit3D::CmsFit3D(int n) : ConstrainedFit(n,6,3+2*n) { }

int CmsFit3D::InitialEstimates(const std::vector<TVectorD> &x,TVectorD &par) {
  par[0] = 0.1;
  par[1] = 0;
  par[2] = 0;
  par[3] = 0;
  par[4] = 0;
  par[5] = 0;
  for(int i=0; i<NTrack(); i++ ){
    double pt = cos(x[i][1])/fabs(x[i][0]);
    double px = pt*cos(x[i][2]);
    double py = pt*sin(x[i][2]);
    double pz = pt*tan(x[i][1]);
    par[2] += x[i][4]/cos(x[i][1]);
    par[3] += px;
    par[4] += py;
    par[5] += pz;
  }
  par[2] /= NTrack();
  return 0;
}

TVectorD CmsFit3D::G(const std::vector<TVectorD> &x,
                      const TVectorD &p) const {
  TVectorD g(NConst());
  double s = p[0];
  double d0 = p[1];
  double z0 = p[2];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  double xv = (s*qx-d0*qy)/qxy;
  double yv = (s*qy+d0*qx)/qxy;
  double zv = z0 + qz*s/qxy;
  g[0] = -qx;
  g[1] = -qy;
  g[2] = -qz;
  for(int i=0; i<NTrack(); i++) {
    double p = 1/fabs(x[i][0]);
    double px = p*cos(x[i][1])*cos(x[i][2]);
    double py = p*cos(x[i][1])*sin(x[i][2]);
    double pz = p*sin(x[i][1]);
    g[0] += px;
    g[1] += py;
    g[2] += pz;
    double s0 = xv*cos(x[i][2]) + yv*sin(x[i][2]);
    g[3+2*i] = x[i][3] + xv*sin(x[i][2]) - yv*cos(x[i][2]);
    g[4+2*i] = x[i][4] - zv*cos(x[i][1]) + s0*sin(x[i][1]);
  }
  //g.Print();
  return g;
}

TMatrixD CmsFit3D::A(const std::vector<TVectorD> &x,
                      const TVectorD &p) const {
  double s = p[0];
  double d0 = p[1];
  double z0 = p[2];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  double xv = (s*qx-d0*qy)/qxy;
  double yv = (s*qy+d0*qx)/qxy;
  double zv = z0 + s*qz/qxy;
  TMatrixD a(5*NTrack(),NConst());
  for(int i=0; i<NTrack(); i++) {
    double pt = cos(x[i][1])/fabs(x[i][0]);
    double px = pt*cos(x[i][2]);
    double py = pt*sin(x[i][2]);
    double pz = pt*tan(x[i][1]);
    a[5*i][0] = -px/x[i][0];
    a[5*i+1][0] = -sin(x[i][1])*cos(x[i][2])/fabs(x[i][0]);
    a[5*i+2][0] = -cos(x[i][1])*sin(x[i][2])/fabs(x[i][0]);
    a[5*i][1] = -py/x[i][0];
    a[5*i+1][1] = -sin(x[i][1])*sin(x[i][2])/fabs(x[i][0]);
    a[5*i+2][1] = cos(x[i][1])*cos(x[i][2])/fabs(x[i][0]);
    a[5*i][2] = -pz/x[i][0];
    a[5*i+1][2] = cos(x[i][1])/fabs(x[i][0]);
    a[5*i+2][3+2*i] = xv*cos(x[i][2])+yv*sin(x[i][2]);
    a[5*i+3][3+2*i] = 1.0;
    a[5*i+1][4+2*i] = zv*sin(x[i][1])+(xv*cos(x[i][2])+yv*sin(x[i][2]))*cos(x[i][1]);
    a[5*i+2][4+2*i] = (-xv*sin(x[i][2])+yv*cos(x[i][2]))*sin(x[i][1]);
    a[5*i+4][4+2*i] = 1.0;
  }
  return a;
}

TMatrixD CmsFit3D::B(const std::vector<TVectorD> &x,
                      const TVectorD &p) const {
  double s = p[0];
  double d0 = p[1];
  double qx = p[3];
  double qy = p[4];
  double qz = p[5];
  double qxy = sqrt(qx*qx+qy*qy);
  double xv = (s*qx-d0*qy)/qxy;
  double yv = (s*qy+d0*qx)/qxy;
  double dxds = qx/qxy;
  double dxdd0 = -qy/qxy;
  double dxdqx = s/qxy-xv*qx/(qxy*qxy);
  double dxdqy = -d0/qxy-xv*qy/(qxy*qxy);
  double dyds = qy/qxy;
  double dydd0 = qx/qxy;
  double dydqx = d0/qxy - yv*qx/(qxy*qxy);
  double dydqy = s/qxy - yv*qy/(qxy*qxy);
  double dzds = qz/qxy;
  double dzdqx = -s*qz*qx/pow(qxy,3);
  double dzdqy = -s*qz*qy/pow(qxy,3);
  double dzdqz = s/qxy;
  double dzdz = 1.0;

  TMatrixD b(NConst(),6);
  b[0][3] = -1.0;
  b[1][4] = -1.0;
  b[2][5] = -1.0;
  for(int i=0; i<NTrack(); i++) {
    b[3+2*i][0] = dxds*sin(x[i][2])-dyds*cos(x[i][2]);
    b[3+2*i][1] = dxdd0*sin(x[i][2])-dydd0*cos(x[i][2]);
    b[3+2*i][3] = dxdqx*sin(x[i][2])-dydqx*cos(x[i][2]);
    b[3+2*i][4] = dxdqy*sin(x[i][2])-dydqy*cos(x[i][2]);
    b[4+2*i][0] = -dzds*cos(x[i][1])+(dxds*cos(x[i][2])+dyds*sin(x[i][2]))*sin(x[i][1]);
    b[4+2*i][1] = (dxdd0*cos(x[i][2])+dydd0*sin(x[i][2]))*sin(x[i][1]);
    b[4+2*i][2] = -dzdz*cos(x[i][1]);
    b[4+2*i][3] = -dzdqx*cos(x[i][1])+(dxdqx*cos(x[i][2])+dydqx*sin(x[i][2]))*sin(x[i][1]);
    b[4+2*i][4] = -dzdqy*cos(x[i][1])+(dxdqy*cos(x[i][2])+dydqy*sin(x[i][2]))*sin(x[i][1]);
    b[4+2*i][5] = -dzdqz*cos(x[i][1]);
  }
  return b;
}

double CmsFit3D::E() const {
  double etot = 0;
  for(int i=0; i<NTrack(); i++) {
    TVectorD x = FitMeasurement(i);
    double p = 1/fabs(x[0]);
    double e = sqrt(p*p+Mass(i)*Mass(i));
    etot += e;
  }
  return etot;
}

double CmsFit3D::M() const {
  return sqrt(E()*E()-Px()*Px()-Py()*Py()-Pz()*Pz());
}

void CmsFit3D::GetVertexMass(unsigned int iv,double &m,double &em) const {
  double q[4] = {0, 0, 0, 0};
  for(int i=0; i<NTrack(); i++) {
    if(Vertex(i) == iv) {
      TVectorD x = FitMeasurement(i);
      double p = 1/fabs(x[0]);
      double px = p*cos(x[1])*cos(x[2]);
      double py = p*cos(x[1])*sin(x[2]);
      double pz = p*sin(x[1]);
      double m = Mass(i);
      double e = sqrt(p*p+m*m);
      q[0] += px;
      q[1] += py;
      q[2] += pz;
      q[3] += e;
    }
  }
  m = sqrt(q[3]*q[3]-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
  // We won't calculate the uncertainty yet...
  em = 0;
}

double CmsFit3D::VertexMass(unsigned int iv) const {
  double m, em;
  GetVertexMass(iv,m,em);
  return m;
}
