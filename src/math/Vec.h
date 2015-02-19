#pragma once
#include "..\sys.h"
#include <iostream>
#include <string>


using namespace std;

template<uint16 dim> 
class Vec
{
public:
	Vec() {
		assert(dim);
		memset(data,0,sizeof(double)*dim); //не очень красиво
	}
	Vec(double first, ...) {
		assert(dim);
		va_list argp;
		va_start(argp, first);
		data[0]=first;
		for (uint16 i=1; i < dim; i++) data[i]=va_arg(argp, double);
		va_end(argp);
	}
	~Vec(void) { }
	double& operator[] (uint16 n);
	const double operator[] (uint16 n) const;
	void display ();
	double lenght ();
	double qlenght ();
	Vec operator+ (const Vec<dim> &op);
	Vec& operator+= (const Vec<dim> &op);
	Vec operator- ();
	Vec& operator= (const Vec<dim> &op);
	Vec operator- (const Vec<dim> &op);
	Vec operator* (const double op);
	string toString ();
	double operator* (const Vec<dim> &op);
	friend std::ostream& operator<< <dim>(std::ostream& stream,const Vec<dim>& obj);
	friend Vec operator* (const double op1, const Vec<dim> &op2);
	double* ptr ();
private:
	double data[dim];
};
//------operator[]------------------------------------------------------
template<uint16 dim>
double& Vec<dim>::operator[] (uint16 n) {
	assert(n < dim);
	return data[n];
}
//-------operator[]-const-----------------------------------------------
template<uint16 dim>
const double Vec<dim>::operator[] (uint16 n) const {
	assert(n < dim);
	return data[n];
}
//--------display----------------------------------------------------
template<uint16 dim>
void Vec<dim>::display () {
	for (uint16 i = 0; i < dim; i++) {
		std::cout << data[i];
		if (i < dim-1) std::cout << ", ";
	}
}
//-------operator+-----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator+(const Vec<dim> &op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]+op.data[i];
	return p;
}
//-------operator+=-----------------------------------------------------
template<uint16 dim>
Vec<dim>& Vec<dim>::operator+= (const Vec<dim> &op) {
	for (uint16 i=0; i < dim; i++) this->data[i]+=op.data[i];
	return *this;
}
//-------operator- ----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator-() {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=-data[i];
	return p;
}
//------operator- -----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator-(const Vec<dim> &op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]-op.data[i];
	return p;
}
//------operator*------------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator*(const double op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]*op;
	return p;
}
//------operator*------------------------------------------------------
template<uint16 dim>
double Vec<dim>::operator*(const Vec<dim> &op) {
	double p=0;
	for (uint16 i=0; i < dim; i++) p+=data[i]*op.data[i];
	return p;
}
//---------qlenght----------------------------------------------------
template<uint16 dim>
double Vec<dim>::qlenght() {
	double p=0;
	for (uint16 i=0; i < dim; i++) p+=data[i]*data[i];
	return p;
}
//---------lenght----------------------------------------------------
template<uint16 dim>
double Vec<dim>::lenght() {
	return sqrt(qlenght());
}
//---------operator*----------------------------------------------------
template<uint16 dim> Vec<dim> operator*(const double op1, const Vec<dim> &op2) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=op2.data[i]*op1;
	return p;
}
//--------operator=------------------------------------------------------
template<uint16 dim> Vec<dim>& Vec<dim>::operator= (const Vec<dim> &op)
{
	memcpy(this->data, op.data, sizeof(double)*dim);
	return *this;
}
//---------operator<<----------------------------------------------------
template<uint16 dim> std::ostream &operator<<(std::ostream &stream, const Vec<dim> &obj) {
	for (uint16 i = 0; i < dim; i++) {
		stream << obj.data[i];
		if (i < dim-1) stream << " "; //TODO: delimeter?
	}
	return stream;
}
//-------------------------------------------------------------
template<uint16 dim> string Vec<dim>::toString ()
{
	string p;
	char buff[100];
	for (uint16 i = 0; i < dim; i++) {
		sprintf_s(buff,100,"%8.5e",this->data[i]);
		p+=buff;
		if (i < dim-1) p+= ", ";
	}
	return p;
}
//-------------------------------------------------------
template<uint16 dim> double* Vec<dim>::ptr ()
{
	return data;
}