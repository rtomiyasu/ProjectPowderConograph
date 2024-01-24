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
#ifndef _RWPARAM_VOID_h_
#define _RWPARAM_VOID_h_
// read_write_func.hh

#include<sstream>

#include "../RietveldAnalysisTypes.hh"
#include "../zerror_type/ZErrorType.hh"
#include "../zerror_type/error_out.hh"
#include "RWParam.hh"


class RWParam_void
{
private:
	RWParamProperty prop;
	void* rwparam;
	const void* rwparam_data;
	
//	template<class T>
//	ZErrorMessage setValue(istream& is, T& ans, const map<string, T>& replace_string_by_value, ostringstream& oss);

	template<class T>
	inline ZErrorMessage setValue(istream& is, ostream* oss);
	template<class T>
	ZErrorMessage setArray(istream& is, ostream* oss);
	template<class T>
	inline ZErrorMessage checkValue() const;
	inline ZErrorMessage checkString() const;
	template<class T>
	ZErrorMessage checkArray() const;

public:
	RWParam_void(const string& label)
		: prop(VOIDDATA, label), rwparam(NULL), rwparam_data(&putVoidTemplate()){};

	template<class T>
	RWParam_void(const pair<RWParamProperty, RWParamData<T> >& arg, T* arg2=NULL)
		: prop(arg.first), rwparam(arg2), rwparam_data(&arg.second){};

	template<class T>
	RWParam_void(const pair<RWParamProperty, RWParamData<T> >& arg, const pair<string, string>& attribute, T* arg2=NULL)
		: prop(arg.first, attribute), rwparam(arg2), rwparam_data(&arg.second){};

	inline const eRWParamType& putType() const { return prop.putType(); };
	inline const string& putLabel() const { return prop.putLabel(); };
	inline string putLongLabel() const { return prop.putLongLabel(); };
	inline const Int4& putMaxMultiNumber() const { return prop.putMaxMultiNumber(); };
	inline const Int4& putMinMultiNumber() const { return prop.putMinMultiNumber(); };
	inline const RWParamProperty& putProperty() const { return prop; };

	inline const void* putParam() const { return rwparam; };
	inline const void* putParamData() const { return rwparam_data; };

	inline void setAttribute(const pair<string, string>& arg){ prop.setAttribute(arg); };

	inline ZErrorMessage setData(istream& is, ostream* oss);
	inline ZErrorMessage checkData() const;

	static const RWParamData<void*>& putVoidTemplate();
};



template<class T>
inline ZErrorMessage RWParam_void::setValue(istream& is, ostream* oss)
{
	T* crwparam = static_cast<T*>(rwparam);
	const RWParamData<T>* crwparam_data = static_cast< const RWParamData<T>* >(rwparam_data);
	string str;
	getfirstword(is, str);
if( oss != NULL ) *oss << str << " ";
	istringstream iss2(str);
	return crwparam_data->replace_string_by_value(iss2, *crwparam);
}


template<class T>
ZErrorMessage RWParam_void::setArray(istream& is, ostream* oss)
{
	const RWParamData< vector<T> >* crwparam_data = static_cast< const RWParamData< vector<T> >* >(rwparam_data);
	vector<T>* crwparam = static_cast< vector<T>* >(rwparam);
	crwparam->clear();

	string str2, str;
	ZErrorType zerr = ZErrorNoError;
	while( zerr == ZErrorNoError )
	{
		zerr = getnewline(is, str2);
		if( !is_blank(str2) ) str += " " + str2;
	}
if( oss != NULL ) *oss << str << " ";
	istringstream iss2(str);
	return crwparam_data->replace_string_by_value(iss2, *crwparam);
}


inline ZErrorMessage RWParam_void::setData(istream& is, ostream* oss)
{
	const eRWParamType& etype = prop.putType();
	if( etype == BOOLFLAG )
	{
		return setValue<bool>(is, oss);
	}
	else if( etype == INT4VALUE )
	{
		return setValue<Int4>(is, oss);
	}
	else if( etype == DVALUE )
	{
		return setValue<Double>(is, oss);
	}
	else if( etype == STRVALUE )
	{
		return setValue<string>(is, oss);
	}
	else if( etype == BOOLARRAY )
	{
		return setArray<bool>(is, oss);
	}
	else if( etype == INT4ARRAY )
	{
		return setArray<Int4>(is, oss);
	}
	else if( etype == DARRAY )
	{
		return setArray<Double>(is, oss);
	}
	return ZErrorMessage(ZErrorUndefined, __FILE__, __LINE__, __FUNCTION__);
}


template<class T>
inline ZErrorMessage RWParam_void::checkValue() const
{ 
	const T* crwparam = static_cast<const T*>(rwparam);
	const RWParamData<T>* crwparam_data = static_cast< const RWParamData<T>* >(rwparam_data);

	ZErrorType zerr = crwparam_data->checkValue(*crwparam);
	if( zerr != ZErrorNoError )
	{
		return nerror_out_range(prop.putLongLabel(),
				num2str<T>(crwparam_data->value_range.begin)+" -- "+num2str<T>(crwparam_data->value_range.end),
				__FILE__, __LINE__, __FUNCTION__);
	}
	return ZErrorMessage();
}


inline ZErrorMessage RWParam_void::checkString() const
{
	const string* crwparam = static_cast<const string*>(rwparam);
	const RWParamData<string>* crwparam_data = static_cast< const RWParamData<string>* >(rwparam_data);

	ZErrorType zerr = crwparam_data->checkArraySize(*crwparam);
	if( zerr != ZErrorNoError )
	{
		return ZErrorMessage(ZErrorArgmentSize, "The length of \""+prop.putLongLabel()+"\" is not within "
				+ num2str<Int4>(crwparam_data->min_array_size)+" -- "+num2str<Int4>(crwparam_data->max_array_size),
				__FILE__, __LINE__, __FUNCTION__);
	}
	return ZErrorMessage();
}


template<class T>
inline ZErrorMessage RWParam_void::checkArray() const
{
	const vector<T>* crwparam = static_cast< const vector<T>* >(rwparam);
	const RWParamData< vector<T> >* crwparam_data = static_cast< const RWParamData< vector<T> >* >(rwparam_data);

	ZErrorType zerr = crwparam_data->checkArraySize(*crwparam);
	if( zerr != ZErrorNoError )
	{
		return ZErrorMessage(ZErrorArgmentSize, "The size of the array \""+prop.putLongLabel()+"\" is not within "
				+ num2str<Int4>(crwparam_data->min_array_size)+" -- "+num2str<Int4>(crwparam_data->max_array_size),
				__FILE__, __LINE__, __FUNCTION__);
	}
	
	zerr = crwparam_data->checkValue(*crwparam);
	if( zerr != ZErrorNoError )
	{
		return nerror_out_range(prop.putLongLabel(),
				vec2str<T>(crwparam_data->value_range.begin)+" -- "+vec2str<T>(crwparam_data->value_range.end),
				__FILE__, __LINE__, __FUNCTION__);
	}
	return ZErrorMessage();
}


inline ZErrorMessage RWParam_void::checkData() const
{
	const eRWParamType& etype = prop.putType();
	if( etype == BOOLFLAG )
	{
		return checkValue<bool>();
	}
	else if( etype == INT4VALUE )
	{
		return checkValue<Int4>();
	}
	else if( etype == DVALUE )
	{
		return checkValue<Double>();
	}
	else if( etype == STRVALUE )
	{
		return checkString();
	}
	else if( etype == BOOLARRAY )
	{
		return checkArray<bool>();
	}
	else if( etype == INT4ARRAY )
	{
		return checkArray<Int4>();
	}
	else if( etype == DARRAY )
	{
		return checkArray<Double>();
	}
	else if( etype == VOIDDATA )
	{
		return ZErrorMessage();
	}
	return ZErrorMessage(ZErrorUndefined, __FILE__, __LINE__, __FUNCTION__);
}


#endif
