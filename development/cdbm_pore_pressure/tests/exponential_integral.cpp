/*********************************************************************
   Returns the exponential integral functions E_n(x) and E_i(x).
   C.A. Bertulani        May/15/2000
*********************************************************************/

typedef double Number;

#include <math.h>
#include<iostream.h>

int main(){
	Number expint(int, Number );
	Number ei(Number );
	Number x;
	int j, n;
	j=1;
	cout << "\n\n Enter n and x.\n";
	cout << "\nWanna check? Note that E_0(x) = exp(-x)/x, and E_n(0) = 1/(n-1), \n";
	cout << "For large x: E_i(x) = exp(x)/x \n";
	cout << "For small x: E_i(x) =  gamma (=0.5772...) + ln(x) \n\n";
	cin >>  n >> x;
    for(;;){
		if(j>10){ cout << "\n My patience is over. Stop, please!\n";
			break;
		}
		if(j!=1){
			cout << "\n\n Enter n and x.\n";
			cin >> n >> x;
		}
		cout << "\n E_n(x): " << expint(n,x);
		cout << "\n E_i(x): " << ei(x);
		j=j+1;
	}
	return 0;
}
/*********************************************************************
   Returns the exponential integral function
   E_n(x) = int_1^infinity e^{-x*t}/t^n dt,     for x > 0.
   C.A. Bertulani        May/15/2000
*********************************************************************/
#define EULER 0.5772156649	/* Euler's constant gamma */
#define MAXIT 100		/* Maximum allowed number of iterations. */
#define FPMIN 1.0e-30	/* close to the smallest representable floting-point number. */
#define EPS 6.0e-8		/* Desired relative error, not smaller than the machine precision. */

Number expint(int n, Number x)
{
	int i,ii,nm1;
	Number a,b,c,d,del,fact,h,psi,ans;

	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
	cerr << "\n Bad arguments in expint";
	else {
		if (n == 0) ans=exp(-x)/x;   /* Special case */
		else {
			if (x == 0.0) ans=1.0/nm1;  /* Another special case */

			else {
				if (x > 1.0) {		/* Lentz's algorithm */
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);	/* Denominators cannot be zero */
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans=h*exp(-x);
							return ans;
						}
					}
					cerr << "\n Continued fraction failed in expint";
				} else {
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);	/* Set first term */
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;  /* Compute psi(n) */
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					cerr << "\n series failed in expint";
				}
			}
		}
	}
	return ans;
}
/************************************************************************
   Returns the exponential integral function
   E_i(x) = - int_x^infinity e^{-t}/t dt = int_(-infinity)^x e^{-t}/t dt,     
   for x > 0.
   C.A. Bertulani        May/15/2000
************************************************************************/

Number ei(Number x)
{
	int k;
	Number fact,prev,sum,term;

	if (x <= 0.0) cerr << "\n Bad argument in ei";
	if (x < FPMIN) return log(x)+EULER;	/* Special case: avoid failure of convergence */
	if (x <= -log(EPS)) {				/* test because of underflow.  */ 
		sum=0.0;				/* Use poer series  */
		fact=1.0;
		for (k=1;k<=MAXIT;k++) {
			fact *= x/k;
			term=fact/k;
			sum += term;
			if (term < EPS*sum) break;
		}
		if (k > MAXIT) cerr << "\n Series failed in ei";
		return sum+log(x)+EULER;
	} else {			/* Use asymptotic series. */
		sum=0.0;		/* Start with second term. */
		term=1.0;
		for (k=1;k<=MAXIT;k++) {
			prev=term;
			term *= k/x;
			if (term < EPS) break;
			/* Since final sum is greater than one, term itself approximates the */
			/* relative error. */
			if (term < prev) sum += term;		/* Still converging: add new term. */
			else {
				sum -= prev;		/* Diverging: subtract previous term and exit. */
				break;
			}
		}
		return exp(x)*(1.0+sum)/x;
	}
}
#undef EPS
#undef EULER
#undef MAXIT
#undef FPMIN
/*********************************************************************/

