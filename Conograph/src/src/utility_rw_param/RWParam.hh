/*
 * The MIT License

   Conograph (powder auto-indexing program)

Copyright (c) <2012> <Ryoko Oishi-Tomiyasu, KEK>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 *
 */
#ifndef _RWPARAM_h_
#define _RWPARAM_h_
// read_write_func.hh

#include<assert.h>
#include<map>

#include "../RietveldAnalysisTypes.hh"
#include "../zerror_type/error_mes.hh"
#include "../utility_data_structure/range.hh"

// IVALUE -> Int4,
// IVALUE_MAX -> "MAX" or Int4,
// DVALUE -> Double,
// DVALUE_MAX -> "MAX" or Int4,
// IARRAY -> vector<Int4>
// DARRAY -> vector<Double>
enum eRWParamType{ BOOLFLAG, INT4VALUE, DVALUE, STRVALUE, BOOLARRAY, INT4ARRAY, DARRAY, VOIDDATA };

class RWParamProperty
{
	friend inline bool IsEqualTag(const RWParamProperty& lhs, const RWParamProperty& rhs){ return lhs.parent_prop == rhs.parent_prop && lhs.label == rhs.label; };

private:
	eRWParamType type;
	string label;
	pair<string, string> attribute;
	const RWParamProperty* parent_prop;
	pair<Int4, Int4> multi_flag;

	string putLongLabel(const Int4& num) const { if( num > 10 || parent_prop==NULL ) return putLabel(); else return parent_prop->putLongLabel(num+1)+"."+putLabel(); };
	bool IsAncestor(const RWParamProperty* arg, const Int4& num) const { if( parent_prop == arg ) return true; if( num > 10 || parent_prop == NULL ) return false; else return parent_prop->IsAncestor(arg, num+1); };

public:

	RWParamProperty(const eRWParamType& typ, const string& lab, const RWParamProperty* arg=NULL,
						const Int4 num1=1, const Int4 num2=1)
			: type(typ), label(lab), attribute("", ""), parent_prop(arg), multi_flag(pair<Int4, Int4>(num1, num2)){ assert(num1 >= 1); assert(num1 >= num2); };

	RWParamProperty(const RWParamProperty& rhs, const pair<string, string>& arg)
			: type(rhs.type), label(rhs.label), attribute(arg), parent_prop(rhs.parent_prop), multi_flag(rhs.multi_flag){};

	const eRWParamType& putType() const { return type; };
	const string& putLabel() const { return label; };
	string putLabelWithAttribute() const { if( attribute.first == "" ) return label; else return label + " " + attribute.first + "=\"" + attribute.second + "\""; };
	string putLongLabel() const { return this->putLongLabel(0); };
	bool IsAncestor(const RWParamProperty* arg) const { return this->IsAncestor(arg, 0); };

	void setAttribute(const pair<string, string>& arg){ attribute = arg; };
	const pair<string, string>& putAttribute() const { return attribute; };
	const RWParamProperty* putParentProperty() const { return parent_prop; };

	inline void setParentProperty(const RWParamProperty* arg) { parent_prop = arg; };

	inline const Int4& putMaxMultiNumber() const { return multi_flag.first; };
	inline const Int4& putMinMultiNumber() const { return multi_flag.second; };
};

template<class T>
class RWParamData
{
private:

public:
	T initial_value;
	ZErrorMessage (*replace_string_by_value)(istream&, T&);
	range<T> value_range;
	range<bool(*)(const T&, const T&)> cmp_func;
	Int4 min_array_size;
	Int4 max_array_size;
	
	RWParamData(const T arg, ZErrorMessage(*arg2)(istream&, T&),
					bool(*cmp_begin)(const T&, const T&), const T range_begin, 
					bool(*cmp_end)(const T&, const T&), const T range_end,
					const Int4& min_asize, const Int4& max_asize);
	
	// Check arg using value_range and cmp_func. 
	ZErrorType checkValue(const T& arg) const;
	// Check arg.size() using initial_array_size, min_array_size and max_array_size.
	ZErrorType checkArraySize(const T& arg) const;
};


template<class T>
RWParamData<T>::RWParamData(const T arg,
		ZErrorMessage(*arg2)(istream&, T&),
		bool(*cmp_begin)(const T&, const T&), const T range_begin, 
		bool(*cmp_end)(const T&, const T&), const T range_end,
		const Int4& min_asize, const Int4& max_asize) 
 : initial_value(arg), 	replace_string_by_value(arg2),
 		value_range(range_begin, range_end), cmp_func(cmp_begin, cmp_end),
		min_array_size(min_asize), max_array_size(max_asize)
{
	assert( cmp_func.begin == NULL || cmp_func.begin(initial_value, value_range.begin) );
	assert( cmp_func.end == NULL || cmp_func.end(initial_value, value_range.end) );
}


// Check arg using value_range and cmp_func. 
template<class T>
ZErrorType RWParamData<T>::checkValue(const T& arg) const
{
	if( cmp_func.begin != NULL )
	{
		if( !cmp_func.begin(arg, value_range.begin) ) return ZErrorOutRange;
	}
	if( cmp_func.end != NULL )
	{
		if( !cmp_func.end(arg, value_range.end) ) return ZErrorOutRange;
	}
	return ZErrorNoError;
}

// Check arg.size() using initial_array_size, min_array_size and max_array_size.
template<class T>
ZErrorType RWParamData<T>::checkArraySize(const T& arg) const
{
	const Int4 isize = arg.size();
	if( 0 <= min_array_size && isize < min_array_size ) return ZErrorArgmentSize;
	if( 0 <= max_array_size && max_array_size < isize ) return ZErrorArgmentSize;
	return ZErrorNoError;
}


#endif
