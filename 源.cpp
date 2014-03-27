#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/tools.h>
#include <windows.h>
#include <set>
#include <iostream>
#include <string>
#pragma comment (lib,"NTL.lib")
using namespace std;
using namespace NTL;
const bool PRIME[] = { 0,
0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 0, 1, 0, 0, 0 };
void IQF_MULTI(ZZ_p u1, ZZ_p v1, ZZ_p u2, ZZ_p v2, ZZ_p& u3, ZZ_p& v3, ZZ n, ZZ base){
	ZZ_p::init(n);
	u3 = u1*u2 + to_ZZ_p(base)*v1*v2;
	v3 = u1*v2 + u2*v1;
}
void IQF_SQUARE(ZZ_p u1, ZZ_p v1, ZZ_p& u2, ZZ_p& v2, ZZ n, ZZ base){
	ZZ_p::init(n);
	IQF_MULTI(u1, v1, u1, v1, u2, v2, n, base);
}
void IQF_EXPON(ZZ_p u1, ZZ_p v1, ZZ_p& u2, ZZ_p& v2, ZZ n, ZZ base, ZZ s){
	ZZ_p::init(n);
	u2 = to_ZZ_p(1); v2 = to_ZZ_p(0); ZZ_p w1 = u1; ZZ_p w2 = v1; ZZ t = s;
	while (t != 0){
		if (t % 2 != 0)IQF_MULTI(u2, v2, w1, w2, u2, v2, n, base);
		t = t / 2; IQF_SQUARE(w1, w2, w1, w2, n, base);
	}
}
bool PERFECT_POWER(ZZ n){
	ZZ k, s, x0, temp; s = int(log(n) / log(10)) + 1;	
	for (k = 2; k <= s; k++){		
		x0 = power(to_ZZ("10"), to_long((s - 1) / k + 1));
		while (power(x0, to_long(k)) > n){
			temp = power(x0, to_long(k - 1));
			x0 = x0 - (temp * x0 - n) / (k*temp) - to_long(1);
		}
		if (power(x0, to_long(k)) == n)return true;
	}
	return false;
}
bool CASE_I(ZZ n)/*n¡Ô1(mod 4)*/{
	ZZ a, k, s;
	for (a = 2;; a++)if (Jacobi(a, n) == -1)break;//generate a
	k = NumTwos(n)*(n - 1); s = ceil(2 * log(log(n)));//generate k s
	if (PowerMod(a, (n - 1) / 2, n) != n - 1)return false;//(1)/a
	if (2 * k > log(n))return true;//(1)/b
	if (PERFECT_POWER(n))return false;//(2)
	ZZ temp1, temp2, m;
	if (s - k<0)temp1 = 0; else temp1 = s - k;
	temp1 = power2_ZZ(to_long(temp1));//cardinality of set S; temp1 =2^max(s-k,0)
	std::set<ZZ> S, S1;
	std::set<ZZ>::iterator p, p1;
	m = 1;
	S.insert(conv<ZZ>("1")); S1.insert(conv<ZZ>("1"));
	temp2 = power2_ZZ(to_long(k));//temp2 =2^k
	while (S.size() < temp1/**/){
		while (S1.find(PowerMod(m, temp2, n)) != S1.end())m++;
		if (m>S.size()*temp2 + 1)return false;
		if (GCD(m, n) > 1)return false;
		ZZ temp3 = power(m, to_long(temp2));//temp3 =m^(2^k)
		for (p1 = S1.begin(); p1 != S1.end(); p1++){
			if (GCD(temp3 - *p1, n) > 1)return false;
		}
			S.insert(m); S1.insert(PowerMod(m, temp2, n));
	}
	for (p = S.begin(); p != S.end(); p++){
		ZZ_p::init(n);
		ZZ_pX poly2 = ZZ_pX(to_long(n), to_ZZ_p(m)) + 1;//poly2 =1+m*x^n
		ZZ_pX poly1 = ZZ_pX(1, to_ZZ_p(m)) + 1;//poly1 =1+m*x
		ZZ_pX module = ZZ_pX(to_long(power2_ZZ(to_long(s))), 1) - to_ZZ_p(a);
		if (PowerMod(poly1%module, n, module) != PowerMod(poly2%module, 1, module))return false;
	}
	return true;
}
bool CASE_II(ZZ n)/*n¡Ô3(mod 4)*/{
	ZZ a, k, t, m, temp;
	for (a = 2;; a++)if (Jacobi(a, n) == -1 && Jacobi((1 - a)%n, n) == -1)break;//generate a	
	k = NumTwos(n)*(n - 1); t = ceil(2 * log(log(n))) + 1;//generate k t	
	if (PowerMod(a, (n - 1) / 2, n) != n - 1)return false;
	ZZ_p::init(n);
	ZZ_p p, q;
	IQF_EXPON(to_ZZ_p(1), to_ZZ_p(1), p, q, n, 1 - a, n);
	if (p != to_ZZ_p(1) || q != to_ZZ_p(-1))return false;
	if (2 * k > log(n))return true;
	if (PERFECT_POWER(n))return false;
	if (t - k > 0)temp = t - k; else temp = 0;
	temp = power2_ZZ(to_long(temp));//temp =2^max(t-k,0)
	for (m = 1; m <= temp; m++)if (GCD(m, n) > 1)return false;
	if (t - k - 1 > 0)temp = t - k - 1; else temp = 0;
	temp = power2_ZZ(to_long(temp));//temp =2^max(t-k-1,0)
	for (m = 1; m <= temp; m++){
		ZZ temp2 = power2_ZZ(to_long(t));//temp2 =2^t
		ZZ_p::init(n);
		ZZ_pX poly1 = ZZ_pX(1, to_ZZ_p(m)) + 1;//poly1 =1+m*x
		ZZ_pX poly2 = ZZ_pX(to_long(n), to_ZZ_p(m)) + 1;//poly2 =1+m*x^n
		ZZ_pX module = ZZ_pX(to_long(temp2 *2), 1) - 2 * ZZ_pX(to_long(temp2), 1) + to_ZZ_p(a);//module x^(2^(t+1))+2*x^(2^t)+a
		if (PowerMod(poly1%module, n, module) != PowerMod(poly2%module, 1, module))return false;
	}
	return true;
}
bool AKS_BER_TEST(ZZ n){
	if (n < 100)return PRIME[to_long(n)];
	else if (!IsOdd(n))return false;
	else if (n % 4 == 1)return CASE_I(n);
	else return CASE_II(n);	
}
int main(){
	ZZ n;
	DWORD E_T; DWORD S_T;
	while (cin >> n){
		S_T = GetTickCount();
		cout << AKS_BER_TEST(n) << ' ';
		E_T = GetTickCount();
		cout << (E_T - S_T) << "ms\n\n" << endl;
	}
	return 0;
}
