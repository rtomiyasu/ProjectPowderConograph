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
#ifndef _I_ReadDATA_h_
#define _I_ReadDATA_h_
// read_write_func.hh

#include <limits>
#include "RWParam.hh"
#include "RWParam_void.hh"
#include "../zerror_type/error_out.hh"

class I_ReadData
{
private:
	virtual void setData(const RWParamProperty& parent_prop,
								vector<RWParam_void>& tray) = 0;
protected:
	static const string brace_opener; // = '<'
	static const string brace_closer; // = '>'

	// On output, all characters in <> are set in the second argument.
	static ZErrorMessage readLabelAll(istream& is, string& label);

	// On output, only the name of tag in <> is set in the second argument.
	static ZErrorMessage readLabel(istream& iss, string& label);

	static ZErrorMessage readToLabelEnd(istream& is, const string& label, string& contents);
	static ZErrorMessage readContents(istream& is, const string& label, string& s);

	static ZErrorMessage removeLabel(istream& is, string& str);
	static ZErrorMessage readAttribute(istream& iss2, map<string, string>& attribute_tray);

	virtual ZErrorMessage checkData(const RWParam_void& param) const;
	virtual ZErrorMessage checkIfDataAreSet(const RWParam_void& param,
											const Int4& set_tray) const;

	ZErrorMessage readData(istream& iss,
						vector<RWParam_void>& param_map, 
						vector<Int4>& set_tray,
						const Int4& outwidth,
						const Int4& stage, ostream* os
					);

	virtual Int4 putOutputNumber(const RWParamProperty& parent_prop) const { return MAX_INT(); };

	template<class T>
	static const map<string, T>& MAP_TO_REPLACE_MIN();
	template<class T>
	static const map<string, T>& MAP_TO_REPLACE_MIN2ZERO();
	template<class T>
	static const map<string, T>& MAP_TO_REPLACE_MAX();

	static map<string, Int4> MAP_TO_REPLACE_NUM_THREAD();

	template<class T>
	static ZErrorMessage setValue(
			const map<string, T>& replace_string_by_value, istream& is, T& ans);

	static ZErrorMessage checkIfDataIsSet(const RWParam_void& param,
								 const Int4& set_tray,
								 const string& attribute_values,
								 const Int4& min_number,
								 const Int4& max_number);
public:
	static Int4 MAX_INT(){ static const Int4 ans = numeric_limits<Int4>::max(); return ans; };
	static Double MAX_DP(){ static const Double ans = numeric_limits<Double>::max(); return ans; };

	template<class T> 
	static bool LE(const T& lhs, const T& rhs){ return lhs <= rhs; };
	template<class T> 
	static bool LT(const T& lhs, const T& rhs){ return lhs < rhs; };
	template<class T> 
	static bool GE(const T& lhs, const T& rhs){ return lhs >= rhs; };
	template<class T> 
	static bool GT(const T& lhs, const T& rhs){ return lhs > rhs; };
	template<class T> 
	static bool GTVec(const vector<T>& lhs, const vector<T>& rhs);

	template<class T>
	static ZErrorMessage REPLACE_NONE(istream&, T&);
	template<class T>
	static ZErrorMessage REPLACE_ENUM(istream& is, T& arg)
	{
		Int4 n;
		is >> n;
		if( is.fail() ) return ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__);
		arg = T(n);
	    return ZErrorMessage();
	};
	template<class T>
	static ZErrorMessage REPLACE_VECTOR_NONE(istream&, vector<T>&);
	template<class T>
	static ZErrorMessage REPLACE_MIN(istream&, T&);
	template<class T>
	static ZErrorMessage REPLACE_MIN2ZERO(istream&, T&);
	template<class T>
	static ZErrorMessage REPLACE_MAX(istream&, T&);

	static ZErrorMessage REPLACE_NUM_THREAD(istream&, Int4&);

	virtual ~I_ReadData(){};

	virtual const string& putTagLabel() const = 0;
	virtual void putParamToReadOrWrite(vector<RWParam_void> & tray) const;

	virtual ZErrorMessageReadingFile readFile(const string& filename, const string& file_label);
	virtual ZErrorMessage readStream(istream& ifs, const string& file_label);
	virtual ZErrorMessage readData(const string& s);
};


template<class T> 
bool I_ReadData::GTVec(const vector<T>& lhs, const vector<T>& rhs)
{ 
	typename vector<T>::const_iterator it=lhs.begin();
	typename vector<T>::const_iterator it2=rhs.begin();
	for(; it<lhs.end() && it2<rhs.end(); it++, it2++)
	{
		if( *it <= *it2 ) return false;
	}
	return true;
}

template<class S, class T>
inline map<S, T> map_term(const S& index, const T& coef)
{
	map<S, T> ans;
	ans.insert(typename map<S, T>::value_type( index, coef ) );
	return ans;
}

template<class T>
const map<string, T>& I_ReadData::MAP_TO_REPLACE_MIN()
{ 
	static const map<string, T> replace_min = map_term<string, T>("MIN", -numeric_limits<T>::max() );
	return replace_min;
};

template<class T>
const map<string, T>& I_ReadData::MAP_TO_REPLACE_MIN2ZERO()
{
	static const map<string, T> replace_min = map_term<string, T>("MIN", 0);
	return replace_min;
};

template<class T>
const map<string, T>& I_ReadData::MAP_TO_REPLACE_MAX()
{
	static const map<string, T> replace_max = map_term<string, T>("MAX", numeric_limits<T>::max() );
	return replace_max;
};



template<class T>
ZErrorMessage I_ReadData::setValue(
		const map<string, T>& replace_string_by_value, istream& is, T& ans)
{
	string str;
	is >> str;
	if( is.fail() ) return ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__);

	typename map<string, T>::const_iterator it_map = replace_string_by_value.find(str);
	if( it_map != replace_string_by_value.end() )
	{
		ans = it_map->second;
	}
	else
	{
	    istringstream iss2(str);
		iss2 >> ans;
		if( iss2.fail() ) return ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__);
	}
	return ZErrorMessage();
}


template<class T>
ZErrorMessage I_ReadData::REPLACE_NONE(istream& is, T& arg)
{
	is >> arg;
	if( is.fail() ) return ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__);
	return ZErrorMessage();
}

template<class T>
ZErrorMessage I_ReadData::REPLACE_VECTOR_NONE(istream& is, vector<T>& arg)
{
	arg.clear();
	T t;
	while( true )
	{
		is >> t;
		if( is.fail() ) return ZErrorMessage();
		arg.push_back(t);
	}
	return ZErrorMessage();
};

template<class T>
ZErrorMessage I_ReadData::REPLACE_MIN(istream& is, T& arg)
{
	return I_ReadData::setValue(MAP_TO_REPLACE_MIN<T>(), is, arg);
};

template<class T>
ZErrorMessage I_ReadData::REPLACE_MIN2ZERO(istream& is, T& arg)
{
	return I_ReadData::setValue(MAP_TO_REPLACE_MIN2ZERO<T>(), is, arg);
};

template<class T>
ZErrorMessage I_ReadData::REPLACE_MAX(istream& is, T& arg)
{
	return I_ReadData::setValue(MAP_TO_REPLACE_MAX<T>(), is, arg);
};


#endif
