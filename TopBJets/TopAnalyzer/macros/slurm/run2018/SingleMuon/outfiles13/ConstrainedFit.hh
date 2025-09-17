#include "TVectorD.h"

#define VERTEX_MASK 0x0f

#define ABZ (1.5e-3*3.8)

class ConstrainedFit {
  int _ncon, _npar, _ntrk, _order, _maxtrk;
  std::vector<TVectorD> _x, _xx;
  std::vector<TMatrixD> _vv;
  std::vector<double> _mass;
  std::vector<unsigned int> _mask;
  TVectorD _par;
  TMatrixD _cov;
  double _chi2;
  int _status;
  double _abz;
public:
  ConstrainedFit(int,int,int);
  virtual ~ConstrainedFit();
  int _verbose;
  double ABz() const { return _abz; }
  int NTrack() const { return _ntrk; }
  virtual int NPar() const { return _npar; }
  virtual int NConst() const { return _ncon; }
  virtual int InitialEstimates(const std::vector<TVectorD> &,TVectorD &);
  virtual void CalculateDerivatives(const std::vector<TVectorD> &,
                                                      const TVectorD &);
  double Chi2() const { return (_status==0)?_chi2:-1; }
  const TVectorD &FitMeasurement(int i) const { return _x[i]; }
  const TVectorD &FitPar() const { return _par; }
  const TMatrixD &FitCov() const { return _cov; }
  double Mass(int i) const { return _mass[i]; }
  virtual void GetVertexMass(unsigned int,double &,double &) const;
  virtual double VertexMass(unsigned int=0) const;
  unsigned int Mask(int i) const { return _mask[i]; }
  unsigned int Vertex(int i) const { return (VERTEX_MASK&_mask[i]); }
  void AddTrack(const double *p,const double *c,double m=0,unsigned int f=0);
  void AddData(const TVectorD &,const TMatrixD &,double m=0,unsigned int f=0);
  void SetVerbose(int v) { _verbose = v; }
  int CheckA();
  int CheckB();
  int Fit(int niter=2);
  virtual TVectorD G(const std::vector<TVectorD> &,
                                const TVectorD &) const = 0;
  virtual TMatrixD A(const std::vector<TVectorD> &,
                               const TVectorD &) const = 0;
  virtual TMatrixD B(const std::vector<TVectorD> &,
                               const TVectorD &) const = 0;
  TVectorD dGdP(int);
  TVectorD dGdX(int,int);
};

class Fit3D : public ConstrainedFit {
  int mass_constraint;
public:
  Fit3D(int);
  ~Fit3D() { }
  int NConst() const { return 3+2*NTrack()+mass_constraint; }
  int InitialEstimates(const std::vector<TVectorD> &,TVectorD &);
  TVectorD G(const std::vector<TVectorD> &,const TVectorD &) const;
  TMatrixD A(const std::vector<TVectorD> &,const TVectorD &) const;
  TMatrixD B(const std::vector<TVectorD> &,const TVectorD &) const;
};

class CmsFit3D : public ConstrainedFit {
public:
  CmsFit3D(int);
  ~CmsFit3D() { }
  int NConst() const { return 3+2*NTrack(); }
  int InitialEstimates(const std::vector<TVectorD> &,TVectorD &);
  TVectorD G(const std::vector<TVectorD> &,const TVectorD &) const;
  TMatrixD A(const std::vector<TVectorD> &,const TVectorD &) const;
  TMatrixD B(const std::vector<TVectorD> &,const TVectorD &) const;
  double S() const { return FitPar()[0]; }
  double ErrS() const { return sqrt(FitCov()[0][0]); }
  double D0() const { return FitPar()[1]; }
  double ErrD0() const { return sqrt(FitCov()[1][1]); }
  double Z0() const { return FitPar()[2]; }
  double ErrZ0() const { return sqrt(FitCov()[2][2]); }
  double Px() const { return FitPar()[3]; }
  double Py() const { return FitPar()[4]; }
  double Pz() const { return FitPar()[5]; }
  double Pt() const { return sqrt(Px()*Px()+Py()*Py()); }
  void GetVertexMass(unsigned int,double &,double &) const;
  double VertexMass(unsigned int=0) const;
  double E() const;
  double M() const;
  double X() const { return (-D0()*Py()+S()*Px())/Pt(); }
  double Y() const { return (D0()*Px()+S()*Py())/Pt(); }
  double Z() const { return Z0()+S()*Pz()/Pt(); }
};
