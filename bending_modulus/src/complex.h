// ver 0.0.1		date 000709

/*
Complex class
–{‘•t‘®‚Ì complex.h “à‚Å’è‹`‚³‚ê‚é•¡‘f”ƒNƒ‰ƒX Complex ‚ÍA
property‚ÉdoubleŒ^‚Ì real ‚Æ imag ‚Ì2¬•ª‚ğ‚à‚¿‘½‚­‚Ì‰‰Zq‚Ì‘½d’è‹`‚ÆŠÖ”‚ğ‚ÂB
ˆÈ‰º‚É‰‰Z‚Ìí—Ş‚Ìˆê——‚ğÚ‚¹‚éB 

‘ã“ü‰‰Zq(’PƒA•¡‡) 

Complex =  Complex
Complex += Complex
Complex -= Complex
Complex *= double
Complex *= Complex
Complex /= double
Complex /= Complex


Zp‰‰Zq(’P€A2€) 

        - Complex
        + Complex
Complex + Complex
Complex - Complex
Complex * double
double  * Complex
Complex * Complex
Complex / double
double  / Complex
Complex / Complex


ŠÖ” 

    // return complex conjugate
Complex conj( Complex& )
    // return absolute value
double  abs( Complex& )
    // return square of absolute value
double  abs2( Complex& )
    // return argument value
double  arg( Complex& )
    // return complex exponential
Complex exp( Complex& )
*/
//////////////////////////////////////////////////////////////////////
#ifndef __COMPLEX_H_INCLUDE
#define __COMPLEX_H_INCLUDE

#include <stdio.h>
#include <math.h>

#define USE_TMPCOMP

//----  •¡‘f”ƒNƒ‰ƒX‚ÌéŒ¾

#ifdef USE_TMPCOMP
class TmpComp
{
public:
  double real, imag;

  TmpComp( void );
  TmpComp( double _real, double _imag=0 );
};
#endif

class Complex
{
public:
  double real;
  double imag;

  Complex( void );
  Complex( double _real, double _imag=0 );

  Complex& operator =( const double& );
  Complex& operator =( const Complex& );
  Complex& operator+=( const Complex& );
  Complex& operator-=( const Complex& );
  Complex& operator*=(       double   );
  Complex& operator*=( const Complex& );
  Complex& operator/=(       double   );
  Complex& operator/=( const Complex& );

#ifdef USE_TMPCOMP
  Complex( const TmpComp& );
  Complex& operator =( const TmpComp& );
  Complex& operator+=( const TmpComp& );
  Complex& operator-=( const TmpComp& );
  Complex& operator*=( const TmpComp& );
  Complex& operator/=( const TmpComp& );
#endif
};

#ifndef USE_TMPCOMP
typedef Complex TmpComp;
#endif


//----  •¡‘f”ƒNƒ‰ƒX‚Ìƒƒ“ƒoŠÖ”

//----  ƒRƒ“ƒXƒgƒ‰ƒNƒ^
inline Complex::Complex( void ){}
inline Complex::Complex( double _real, double _imag )
{
  real = _real;    imag = _imag;
}

#ifdef USE_TMPCOMP
inline TmpComp::TmpComp( void ){}
inline TmpComp::TmpComp( double _real, double _imag )
{
  real = _real;    imag = _imag;
}
inline Complex::Complex( const TmpComp& t )
{
  real = t.real;    imag = t.imag;
}
#endif

//----  ‘ã“ü‚Ì’è‹`
inline Complex&  Complex::operator = ( const double& c )
{
  real = c;   imag = 0.0;
  return *this;
}
inline Complex&  Complex::operator = ( const Complex& c )
{
  real = c.real;   imag = c.imag;
  return *this;
}
#ifdef USE_TMPCOMP
inline Complex&  Complex::operator = ( const TmpComp& t )
{
  real = t.real;   imag = t.imag;
  return *this;
}
#endif

//----  ’P€‰‰Zq‚Ì’è‹`
inline TmpComp  operator - ( const Complex& c )
{
  return( TmpComp( -c.real, -c.imag ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator - (       TmpComp& t )
{
  t.real *= -1;    t.imag *= -1;
  return( t );
}
inline Complex& operator + (       Complex& c )
{
  return( c );
}
inline TmpComp& operator + (       TmpComp& t )
{
  return( t );
}
#endif

//----  •¡‘f” += •¡‘f”‚Ì’è‹`
inline Complex& Complex::operator+=( const Complex& c )
{
  real += c.real;    imag += c.imag;
  return( *this );  
}
#ifdef USE_TMPCOMP
inline Complex& Complex::operator+=( const TmpComp& t )
{
  real += t.real;    imag += t.imag;
  return( *this );  
}
#endif

//----  •¡‘f” -= •¡‘f”‚Ì’è‹`
inline Complex& Complex::operator-=( const Complex& c )
{
  real -= c.real;    imag -= c.imag;
  return( *this );  
}
#ifdef USE_TMPCOMP
inline Complex& Complex::operator-=( const TmpComp& t )
{
  real -= t.real;    imag -= t.imag;
  return( *this );  
}
#endif

//----  •¡‘f” *= À”‚Ì’è‹`
inline Complex& Complex::operator*=( double k )
{
  real *= k;    imag *= k;
  return( *this );
}

//----  •¡‘f” *= •¡‘f”‚Ì’è‹`
inline Complex& Complex::operator*=( const Complex& c )
{
  const double temp = real*c.real - imag*c.imag;
  imag = real*c.imag + imag*c.real;
  real = temp;
  return( *this );
}
#ifdef USE_TMPCOMP
inline Complex& Complex::operator*=( const TmpComp& t )
{
  const double temp = real*t.real - imag*t.imag;
  imag = real*t.imag + imag*t.real;
  real = temp;
  return( *this );
}
#endif

//----  •¡‘f” /= À”‚Ì’è‹`
inline Complex& Complex::operator/=( double k )
{
  real /= k;    imag /= k;
  return( *this );
}

//----  •¡‘f” /= •¡‘f”‚Ì’è‹`
inline Complex& Complex::operator/=( const Complex& c )
{
  const double fact = 1.0/(c.real*c.real + c.imag*c.imag);
  const double temp = (real*c.real + imag*c.imag) * fact;
  imag = (imag*c.real - real*c.imag) * fact;
  real = temp;
  return( *this );
}
#ifdef USE_TMPCOMP
inline Complex& Complex::operator/=( const TmpComp& t )
{
  const double fact = 1.0/(real*t.real + imag*t.imag);
  const double temp = (real*t.real + imag*t.imag) * fact;
  imag = (imag*t.real - real*t.imag) * fact;
  real = temp;
  return( *this );
}
#endif

//----  •¡‘f” + À”‚Ì’è‹`
inline TmpComp  operator+( const Complex& c, const double& real )
{
  return( TmpComp( c.real+real, c.imag )  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator+(       TmpComp&  t, const double& real )
{
  t.real += real;
  return( t );
}
#endif

//----  À” + •¡‘f”‚Ì’è‹`
inline TmpComp  operator+( const double& real, const Complex& c )
{
  return( TmpComp( real+c.real, c.imag )  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator+( const double& real,      TmpComp&  t )
{
  t.real += real;
  return( t );
}
#endif


//----  •¡‘f” + •¡‘f”‚Ì’è‹`
inline TmpComp  operator+( const Complex& c1, const Complex& c2 )
{
  return( TmpComp( c1.real+c2.real, c1.imag+c2.imag)  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator+(       TmpComp&  t, const Complex&  c )
{
  t.real += c.real;      t.imag += c.imag;
  return( t );
}
inline TmpComp& operator+( const Complex&  c,       TmpComp&  t )
{
  t.real += c.real;      t.imag += c.imag;
  return( t );
}
inline TmpComp& operator+(       TmpComp& t1, const TmpComp& t2 )
{
  t1.real += t2.real;    t1.imag += t2.imag;
  return( t1 );
}
#endif

//----  •¡‘f” - À”‚Ì’è‹`
inline TmpComp  operator-( const Complex& c, const double& real )
{
  return( TmpComp( c.real-real, c.imag )  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator-(       TmpComp&  t, const double& real )
{
  t.real -= real;
  return( t );
}
#endif

//----  À” - •¡‘f”‚Ì’è‹`
inline TmpComp  operator-( const double& real, const Complex& c )
{
  return( TmpComp( real-c.real, -c.imag )  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator-( const double& real,      TmpComp&  t )
{
  t.real *= -1;    t.imag *= -1;
  t.real += real;
  return( t );
}
#endif

//----  •¡‘f” - •¡‘f”‚Ì’è‹`
inline TmpComp  operator-( const Complex& c1, const Complex& c2 )
{
  return( TmpComp( c1.real-c2.real, c1.imag-c2.imag)  );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator-(       TmpComp&  t, const Complex&  c )
{
  t.real -= c.real;      t.imag -= c.imag;
  return( t );
}
inline TmpComp& operator-( const Complex&  c,       TmpComp&  t )
{
  t.real *= -1;          t.imag *= -1;
  t.real += c.real;      t.imag += c.imag;
  return( t );
}
inline TmpComp& operator-(       TmpComp& t1, const TmpComp& t2 )
{
  t1.real -= t2.real;    t1.imag -= t2.imag;
  return( t1 );
}
#endif


//----  À” * •¡‘f”‚Ì’è‹`
inline TmpComp  operator*( double k, const Complex& c )
{
  return( TmpComp( k*c.real, k*c.imag ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator*( double k,       TmpComp& t )
{
  t.real *= k;      t.imag *= k;
  return( t );
}
#endif

//----  •¡‘f” * À”‚Ì’è‹`
inline TmpComp  operator*( const Complex& c, double k )
{
  return( TmpComp( c.real*k, c.imag*k ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator*(       TmpComp& t, double k )
{
  t.real *= k;      t.imag *= k;
  return( t );
}
#endif

//----  •¡‘f” / À”‚Ì’è‹`
inline TmpComp  operator/( const Complex& c, double k )
{
  const double inv_k = 1.0/k;
  return( TmpComp( c.real*inv_k, c.imag*inv_k ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator/(       TmpComp& t, double k )
{
  const double inv_k = 1.0/k;
  t.real *= inv_k;      t.imag *= inv_k;
  return( t );
}
#endif

//----  À” / •¡‘f”‚Ì’è‹`
inline TmpComp  operator/( double k, const Complex& c )
{
  const double fact   = 1.0/(c.real*c.real + c.imag*c.imag);
  return( TmpComp( k*c.real*fact, -k*c.imag*fact ) );
}
#ifdef USE_TMPCOMP
inline TmpComp  operator/( double k,       TmpComp& t )
{
  const double fact   = 1.0/(t.real*t.real + t.imag*t.imag);
  t.real =  k*t.real*fact;
  t.imag = -k*t.imag*fact;
  return( t );
}
#endif

//----  •¡‘f” * •¡‘f”‚Ì’è‹`
inline TmpComp  operator*( const Complex& c1, const Complex& c2 )
{
  return( TmpComp( c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator*(       TmpComp&  t, const Complex&  c )
{
  const double t_real = t.real;
  t.real = t_real*c.real - t.imag*c.imag;
  t.imag = t_real*c.imag + t.imag*c.real;
  return( t );
}
inline TmpComp& operator*( const Complex&  c,       TmpComp&  t )
{
  const double t_real = t.real;
  t.real = t_real*c.real - t.imag*c.imag;
  t.imag = t_real*c.imag + t.imag*c.real;
  return( t );
}
inline TmpComp& operator*(       TmpComp& t1, const TmpComp& t2 )
{
  const double t1_real = t1.real;
  t1.real = t1_real*t2.real - t1.imag*t2.imag;
  t1.imag = t1_real*t2.imag + t1.imag*t2.real;
  return( t1 );
}
#endif

//----  •¡‘f” / •¡‘f”‚Ì’è‹`
inline TmpComp  operator/( const Complex& c1, const Complex& c2 )
{
  const double fact   = 1.0/(c2.real*c2.real + c2.imag*c2.imag);
  return( TmpComp( (c1.real*c2.real + c1.imag*c2.imag) * fact,
                   (c1.imag*c2.real - c1.real*c2.imag) * fact ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& operator/(       TmpComp&  t, const Complex&  c )
{
  const double t_real = t.real;
  const double fact   = 1.0/(c.real*c.real + c.imag*c.imag);
  t.real = (t_real*c.real + t.imag*c.imag) * fact;
  t.imag = (t.imag*c.real - t_real*c.imag) * fact;
  return( t );
}
inline TmpComp& operator/( const Complex&  c,       TmpComp&  t )
{
  const double t_real = t.real;
  const double fact   = 1.0/(t.real*t_real + t.imag*t.imag);
  t.real = (c.real*t_real + c.imag*t.imag) * fact;
  t.imag = (c.imag*t_real - c.real*t.imag) * fact;
  return( t );
}
inline TmpComp& operator/(       TmpComp& t1, const TmpComp& t2 )
{
  const double t1_real= t1.real;
  const double fact   = 1.0/(t2.real*t2.real + t2.imag*t2.imag);
  t1.real = (t1_real*t2.real + t1.imag*t2.imag) * fact;
  t1.imag = (t1.imag*t2.real - t1_real*t2.imag) * fact;
  return( t1 );
}
#endif

inline double Real( const Complex& c )
{
  return( c.real );
}
inline double Real( const TmpComp& c )
{
  return( c.real );
}
inline double Imag( const Complex& c )
{
  return( c.imag );
}
inline double Imag( const TmpComp& c )
{
  return( c.imag );
}


//----  •¡‘f‹¤–ğ‚Ì’è‹`
inline TmpComp conj( const Complex& c )
{
  return( TmpComp( c.real, -c.imag ) );
}

//----  â‘Î’l‚Ì’è‹`
inline double abs( const Complex& c )
{
#ifdef WIN32
	return( _hypot( c.real, c.imag ) );
#else
	return( hypot( c.real, c.imag ) );
#endif
}

//----  â‘Î’l‚Ì2æ‚Ì’è‹`
inline double norm( const Complex& c )
{
    return( c.real*c.real+c.imag*c.imag );
}
inline double abs2( const Complex& c )
{
    return( c.real*c.real+c.imag*c.imag );
}

//----  •ÎŠp‚Ì’è‹`
inline double arg( const Complex& c )
{
  return( atan2( c.imag, c.real ) );
}

//----  w”ŠÖ”‚Ì’è‹`
inline TmpComp  exp( const Complex& c )
{
  const double R = exp(c.real);
  return( TmpComp( R*cos(c.imag), R*sin(c.imag) ) );
}
#ifdef USE_TMPCOMP
inline TmpComp& exp( TmpComp& t )
{
  const double R = exp(t.real);
  t.real = R*cos(t.imag);
  t.imag = R*sin(t.imag);
  return( t );
}
#endif

//----  w”ŠÖ”‚Ì’è‹` expCI(r) == exp(CI*r)
inline TmpComp  expCI( const double& d )
{
  return( TmpComp( cos(d), sin(d) ) );
}

//----  2æ‚Ì’è‹`
template <class T>
inline T sqr( T c )
{
  return c*c;
}
#ifdef USE_TMPCOMP
inline TmpComp& sqr(       TmpComp& t )
{
  const double t_real = t.real;
  const double t_imag = t.imag;
  t.real = t_real*t_real - t_imag*t_imag;
  t.imag = t_real*t_imag + t_imag*t_real;
  return( t );
}
#endif



//----  •¡‘f”‚Ì¬•ª‚Ì•\¦
inline void print( const Complex& c )
{
  printf("(%+f %+f)\n", c.real, c.imag );
}


//---- À”‚Ì’PˆÊ
const Complex Cr1(1.0, 0.0);

//---- ‹•”‚Ì’PˆÊ
const Complex Ci1(0.0, 1.0);

//---- ‹•”‚Ì0
const Complex C0(0.0, 0.0);


#endif  // __COMPLEX_H_INCLUDE

