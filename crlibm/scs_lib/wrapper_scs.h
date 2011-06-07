extern "C"
{
#include "scs.h"
}
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream.h>

#ifndef SCS_WRAPPER_CPP
#define SCS_WRAPPER_CPP 1

extern "C" int scs_cmp_mant(scs_ptr , scs_ptr);

class Scs{
 private:
    scs scsnb; 

 public:
    /* Constructors */
    Scs(){}
    Scs(double &d){scs_set_d(&scsnb, d);}
    Scs(int &i){scs_set_si(&scsnb, i);}
    Scs(const Scs& nb);
    ~Scs(){}; 

    /* Mutators */
    void scs_set_sign (int i)     {this->scsnb.sign = i;}
    void scs_set_index(int i)     {this->scsnb.index = i;}
    void scs_set_excep(double d)  {this->scsnb.exception.d = d;}
    void scs_set_words(unsigned int word[SCS_NB_WORDS]){for(int i=0; i<SCS_NB_WORDS; i++) this->scsnb.h_word[i] = word[i];}

    /* Cast */
    operator double();
    operator int();


    /* Negation */
    Scs operator-() const;
    friend Scs fabs(const Scs &a);

    /* Assignment */
    Scs &operator=(const Scs &nb1);
    Scs &operator=(const double nb1);
    Scs &operator=(const int i){scs_set_si(&scsnb, i);}

    /* Addition */
    friend Scs operator+(Scs &nb1,Scs &nb2);
    friend Scs operator+(Scs &nb1, const double &nb2);
    friend Scs operator+(const double &nb1,Scs &nb2);
    void operator+=(Scs &nb);
    void operator+=(const double nb);

    /* Subtraction */
    friend Scs operator-(Scs &nb1,Scs &nb2);
    friend Scs operator-(Scs &nb1, const double &nb2);
    friend Scs operator-(const double &nb1,Scs &nb2);
    void operator-=(Scs &nb);
    void operator-=(const double nb);

    /* Multiplication */
    friend Scs operator*(Scs &nb1,Scs &nb2);
    friend Scs operator*(Scs &nb1, const double &nb2);
    friend Scs operator*(const double &nb1,Scs &nb2);
    friend Scs operator*(Scs &nb1, const int &nb2);
    friend Scs operator*(const int &nb1,Scs &nb2);
    void operator*=(Scs &nb);
    void operator*=(const double nb);
    void operator*=(const int nb);

    /* Multiplication */
    friend Scs operator/(Scs &nb1,Scs &nb2);
    friend Scs operator/(Scs &nb1, const double &nb2);
    friend Scs operator/(const double &nb1,Scs &nb2);
    void operator/=(Scs &nb);
    void operator/=(const double nb);

    /* Comparisons */
    friend bool operator==(Scs &nb1,Scs &nb2);
    friend bool operator!=(Scs &nb1,Scs &nb2);
    friend bool operator<=(Scs &nb1,Scs &nb2);
    friend bool operator>=(Scs &nb1,Scs &nb2);
    friend bool operator<(Scs &nb1,Scs &nb2);
    friend bool operator>(Scs &nb1,Scs &nb2);

    /* Random Number */
    Scs rand(void);

    /* Input/Output */
    friend ostream& operator<<(ostream &s, const Scs &a);
    friend istream& operator>>(istream &s, Scs &a);

};




/**************
 * CONSTRUCTOR
 **************/
Scs::Scs(const Scs& nb){
  unsigned int i;
  
  for(i=0; i<SCS_NB_WORDS; i++)
    this->scsnb.h_word[i] = nb.scsnb.h_word[i];

  this->scsnb.exception.d = nb.scsnb.exception.d;
  this->scsnb.index = nb.scsnb.index;
  this->scsnb.sign = nb.scsnb.sign;
}

/**************
 * CAST
 **************/
inline Scs::operator double() {
    double d;
    scs_get_d(&d, &(this->scsnb));
    return d;
}
inline Scs::operator int() {
    double d;
    scs_get_d(&d, &(this->scsnb));
    return ((int)d);
}

 

/**************    
 * ASSIGNATION
 **************/
inline Scs &Scs::operator=(const Scs& nb){
  unsigned int i;
  
  for(i=0; i<SCS_NB_WORDS; i++)
    scsnb.h_word[i] = nb.scsnb.h_word[i];

  scsnb.exception.d = nb.scsnb.exception.d;
  scsnb.index = nb.scsnb.index;
  scsnb.sign = nb.scsnb.sign;
  
  return *this;
}
inline Scs &Scs::operator=(const double nb){
    scs_set_d(&(this->scsnb), nb);
    return *this;
}
inline Scs fabs(const Scs &a){
    Scs res(a);
    res.scsnb.sign = 1;
    return res;
}


/************
 * ADDITION
 ************/
inline Scs operator+(Scs &nb1,Scs &nb2){
    Scs res;    scs_add(&(res.scsnb), &(nb1.scsnb), &(nb2.scsnb));
    return res;
}
inline Scs operator+(Scs &nb1,const double &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb2);
    scs_add(&(res.scsnb), &(nb1.scsnb), &(op.scsnb));
    return res;
}
inline Scs operator+(const double &nb1, Scs &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb1);
    scs_add(&(res.scsnb), &(nb2.scsnb), &(op.scsnb));
    return res;
}
void inline Scs::operator+=(Scs &nb) {
    scs_add(&(this->scsnb), &(this->scsnb), &(nb.scsnb));
}
void inline Scs::operator+=(const double nb) {
    Scs op;
    scs_set_d(&(op.scsnb), nb);
    scs_add(&(this->scsnb), &(this->scsnb), &(op.scsnb));
}




/**************
 * SUBTRACTION
 **************/
inline Scs operator-(Scs &nb1,Scs &nb2){
    Scs res;    scs_sub(&(res.scsnb), &(nb1.scsnb), &(nb2.scsnb));
    return res;
}
inline Scs operator-(Scs &nb1,const double &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb2);
    scs_sub(&(res.scsnb), &(nb1.scsnb), &(op.scsnb));
    return res;
}
inline Scs operator-(const double &nb1, Scs &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb1);
    scs_sub(&(res.scsnb), &(nb2.scsnb), &(op.scsnb));
    return res;
}
void inline Scs::operator-=(Scs &nb) {
    scs_sub(&(this->scsnb), &(this->scsnb), &(nb.scsnb));
}
void inline Scs::operator-=(const double nb) {
    Scs op;
    scs_set_d(&(op.scsnb), nb);
    scs_sub(&(this->scsnb), &(this->scsnb), &(op.scsnb));
}



/*****************
 * MULTIPLICATION
 *****************/
inline Scs operator*(Scs &nb1,Scs &nb2){
    Scs res;    scs_mul(&(res.scsnb), &(nb1.scsnb), &(nb2.scsnb));
    return res;
}
inline Scs operator*(Scs &nb1,const double &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb2);
    scs_mul(&(res.scsnb), &(nb1.scsnb), &(op.scsnb));
    return res;
}
inline Scs operator*(const double &nb1, Scs &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb1);
    scs_mul(&(res.scsnb), &(nb2.scsnb), &(op.scsnb));
    return res;
}
inline Scs operator*(Scs &nb1, const int &nb2){
    Scs res;
    scs_set(&(res.scsnb), &(nb1.scsnb));
    scs_mul_ui(&(res.scsnb), nb2);
    return res;
}
inline Scs operator*(const int &nb1, Scs &nb2){
    Scs res;
    scs_set(&(res.scsnb), &(nb2.scsnb));
    scs_mul_ui(&(res.scsnb), nb1);
    return res;
}
void inline Scs::operator*=(Scs &nb) {
    scs_mul(&(this->scsnb), &(this->scsnb), &(nb.scsnb));
}
void inline Scs::operator*=(const double nb) {
    Scs op;
    scs_set_d(&(op.scsnb), nb);
    scs_mul(&(this->scsnb), &(this->scsnb), &(op.scsnb));
}
void inline Scs::operator*=(const int nb) {
    scs_mul_ui(&(this->scsnb), nb);
}



/*****************
 * DIVISION
 *****************/
inline Scs operator/(Scs &nb1,Scs &nb2){
    Scs res;    scs_div(&(res.scsnb), &(nb1.scsnb), &(nb2.scsnb));
    return res;
}
inline Scs operator/(Scs &nb1, const double &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb2);
    scs_div(&(res.scsnb), &(nb1.scsnb), &(op.scsnb));
    return res;
}
inline Scs operator/(const double &nb1, Scs &nb2){
    Scs res, op;
    scs_set_d(&(op.scsnb), nb1);
    scs_div(&(res.scsnb), &(nb2.scsnb), &(op.scsnb));
    return res;
}
void inline Scs::operator/=(Scs &nb) {
    scs_div(&(this->scsnb), &(this->scsnb), &(nb.scsnb));
}
void inline Scs::operator/=(const double nb) {
    Scs op;
    scs_set_d(&(op.scsnb), nb);
    scs_div(&(this->scsnb), &(this->scsnb), &(op.scsnb));
}



/*************
 * COMPARISON
 *************/
inline bool operator==(Scs &nb1, Scs &nb2){
  unsigned int i;
  bool b=1;
  
  for(i=0; i<SCS_NB_WORDS; i++) 
      if (nb1.scsnb.h_word[i] == nb2.scsnb.h_word[i]) b=0;
  
  return ((nb1.scsnb.exception.d == nb2.scsnb.exception.d)&&
	  (nb1.scsnb.index == nb2.scsnb.index)&&
	  (nb1.scsnb.sign  == nb2.scsnb.sign)&& b);
}
inline bool operator!=(Scs &nb1, Scs &nb2){
    return !(nb1==nb2);
}
inline bool operator<=(Scs &nb1, Scs &nb2){
    return ((nb1.scsnb.exception.d <= nb2.scsnb.exception.d)&&
	    (nb1.scsnb.sign  <= nb2.scsnb.sign)&&
	    ((nb1.scsnb.index < nb2.scsnb.index)||
	     ((nb1.scsnb.index == nb2.scsnb.index)&&
	      (scs_cmp_mant(&(nb1.scsnb), &(nb2.scsnb))<=0))));
}
inline bool operator>=(Scs &nb1, Scs &nb2){
    return ((nb1.scsnb.exception.d >= nb2.scsnb.exception.d)&&
	    (nb1.scsnb.sign  >= nb2.scsnb.sign)&&
	    ((nb1.scsnb.index > nb2.scsnb.index)||
	     ((nb1.scsnb.index == nb2.scsnb.index)&&
	      (scs_cmp_mant(&(nb1.scsnb), &(nb2.scsnb))>=0))));
}
inline bool operator<(Scs &nb1, Scs &nb2){
    return ((nb1.scsnb.exception.d <= nb2.scsnb.exception.d)&&
	    (nb1.scsnb.sign  <= nb2.scsnb.sign)&&
	    ((nb1.scsnb.index < nb2.scsnb.index)||
	     ((nb1.scsnb.index == nb2.scsnb.index)&&
	      (scs_cmp_mant(&(nb1.scsnb), &(nb2.scsnb))<0))));
}
inline bool operator>(Scs &nb1, Scs &nb2){
    return ((nb1.scsnb.exception.d >= nb2.scsnb.exception.d)&&
	    (nb1.scsnb.sign  >= nb2.scsnb.sign)&&
	    ((nb1.scsnb.index > nb2.scsnb.index)||
	     ((nb1.scsnb.index == nb2.scsnb.index)&&
	      (scs_cmp_mant(&(nb1.scsnb), &(nb2.scsnb))>0))));
}




/****************
 * RANDOM NUMBER
 ****************/
inline Scs Scs::rand(void){
    scs_rand(&(this->scsnb), 200); 
    return *this;
}




/*************** 
 * OUTPUT (in hexadecimal)
 ***************/ 
ostream &operator<<(ostream &os, const Scs &a){
    Scs aa, p, zer;
    double d;
    char buffer[10];
    int e, exposant;
    bool bb;

    if (a.scsnb.exception.d != 1.){
	os << (double)a.scsnb.exception.d;
    }else {
	if (a.scsnb.sign == -1)
	    os << '-';
	
	aa = fabs(a);

	/* Compute the exposant in radix 16 */
	d = ((a.scsnb.index)*SCS_NB_BITS)/4;
	e = 4*(int)floor(d);

	p = 1;
	p.scsnb.index = (int)floor(((double)e)/SCS_NB_BITS);
	p.scsnb.h_word[0] = 1 << e - p.scsnb.index*SCS_NB_BITS;
	exposant = (int)floor(d);
	p /= 16;
	exposant--;
	while(p <= aa){
	    p *= 16;
	    exposant++;
	}
	p /= 16;
	exposant--;
	
	/* Extract digits */
	aa = aa / p;
	sprintf(buffer,"%x", aa.scsnb.h_word[0]);
	os << buffer << ".";
	aa.scsnb.h_word[0] = 0;
	aa *= 16;
	
	bb = 1;
	while(bb){
	    sprintf(buffer,"%x", aa.scsnb.h_word[0]);
	    os << buffer;
	    aa.scsnb.h_word[0] = 0;
	    aa *= 16;

	    bb = 0;
	    for(int i=0; i<SCS_NB_WORDS; i++)
		if (aa.scsnb.h_word[i] != 0) bb=1;
	}

	/* Write the exponent */
	os << " x16^(" << exposant <<")";
    } 
    
    return os;
}


/*************** 
 * INPUT (in decimal)
 ***************/ 
istream& operator>>(istream &is, Scs &a){
    char c;
    int nd = 0;
    int point = -1;
    int ex;
    bool done = false;
    Scs r;

    r = 0;

    /* Skip any leading spaces */
    do{
	is>>c;
    }while (c == ' ');

    /* Read sign, digits, and exponent */
    while (!done && (c != '\0')) {
	if (c >= '0' && c <= '9') {
	    int d = c - '0';
	    r *= 10.0;
	    r += d;
	    nd++;
	} else {
	    switch (c) {
	    case '.':
		point = nd;
		break;
	    case '-':
	    case '+':
		if (nd > 0){
		    a = 0;
		    done = true;
		    point = -1;
		    ex = 0;
		}
		a.scsnb.sign = (c == '-') ? -1 : 1;
		break;
	    case 'E':
	    case 'e':
		is >> ex;
		done = true;
		break;
	    default:
		a = 0;
		done = true;
		point = -1;
		ex = 0;
	    }
	}
	is>>c;
    }

    if (point >= 0) 
	ex -= (nd - point);
    
    
    if (ex != 0) {
	if (ex > 0)
	    for(int i=0; i<ex; i++)
		r *= 10;
	if (ex < 0){
	    Scs inv_ten, ten;
	    ten = 10;
	    scs_inv(&(inv_ten.scsnb), &(ten.scsnb) );
	    for(int i=0; i>ex; i--)
		r *= inv_ten;
	}
    }
}
#endif


