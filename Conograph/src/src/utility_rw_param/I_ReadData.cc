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
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <fstream>
#include "I_ReadData.hh"
#include "../zlog/zlog.hh"

//const string I_ReadData::comment_opener = "<!--";
//const string I_ReadData::comment_closer = "-->";
const string I_ReadData::brace_opener = "<";
const string I_ReadData::brace_closer = ">";

ZErrorMessage I_ReadData::readLabelAll(istream& is, string& label)
{
	ZErrorMessage zerr;
	zerr = getdelim(is, label, brace_opener );
	if( zerr.putErrorType() != ZErrorNoError )
	{
    	return zerr;
	}

	if( is.peek() == '!')
	{
		char c;
		is.get(c);

		if( is.peek() == '[')	// ![CDATA[
		{
			zerr = getdelim(is, label, "]]" + brace_closer );
		}
		else // if( is.peek() == '-')	// !--
		{
			zerr = getdelim(is, label, "--" + brace_closer );
		}
		label = "!"+ label;
		return zerr;
	}
	else if( is.peek() == '?')
	{
		return getdelim(is, label, "?" + brace_closer );
	}

	return getdelim(is, label, brace_closer );
}


// Read functions.
ZErrorMessage I_ReadData::readLabel(istream& is, string& label)
{
	ZErrorMessage zerr = readLabelAll(is, label);
	if( zerr.putErrorType() != ZErrorNoError )
	{
    	return zerr;
	}

	istringstream iss(label);
	getfirstword(iss, label);
	return ZErrorMessage();
}



// Read functions.
ZErrorMessage I_ReadData::readToLabelEnd(istream& is, const string& label, string& s)
{
	const string end_label = brace_opener + "/"+label + brace_closer;
	ZErrorMessage zerr = getdelim(is, s,  end_label);
    if( zerr.putErrorType() != ZErrorNoError )
    {
    	return zerr;
    }

    istringstream iss2(s);
    string str;
    Int4 count = 0;
	zerr = getdelim(iss2, str,  brace_opener + label + brace_closer);
    while ( zerr.putErrorType() == ZErrorNoError )
    {
    	count++;
    	zerr = getdelim(iss2, str,  brace_opener + label + brace_closer);
    }

    while( count > 0 )
    {
    	zerr = getdelim(is, str,  end_label);
        if( zerr.putErrorType() != ZErrorNoError )
        {
        	return zerr;
        }
        s += end_label + str;
    	count--;
    }
	return ZErrorMessage();
}



ZErrorMessage I_ReadData::removeLabel(istream& is, string& str)
{
	str.clear();
	string str2;
	ZErrorMessage zerr;
	while( true )
	{
		zerr = getdelim(is, str2, brace_opener );
		str += str2;
		if( zerr.putErrorType() != ZErrorNoError ) break;

		if( is.peek() == '!')
		{
			zerr = getdelim(is, str2, "--" + brace_closer );
		}
		else if( is.peek() == '?')
		{
			zerr = getdelim(is, str2, "?" + brace_closer );
		}
		else zerr = getdelim(is, str2, brace_closer );
		if( zerr.putErrorType() != ZErrorNoError )
		{
	    	return zerr;
		}
		str += " ";
	}
	return ZErrorMessage();
}


ZErrorMessage I_ReadData::readAttribute(istream& iss2, map<string, string>& attribute_tray)
{
	attribute_tray.clear();

	ZErrorMessage zerr;
	string lhs, rhs, str2;
	zerr = getdelim(iss2, lhs, "=");
	while( zerr.putErrorType() == ZErrorNoError )
	{
		zerr = getdelim(iss2, str2, "\"");
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;
		zerr = getdelim(iss2, rhs, "\"");
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;

		istringstream iss3(lhs);
		getfirstword(iss3, lhs);
		if( lhs.empty() ) return ZErrorMessage(ZErrorFileFormatBroken, "No word is written before '='", __FILE__, __LINE__, __FUNCTION__);

		const map<string, string>::iterator it=attribute_tray.find(lhs);
		if( it != attribute_tray.end() )
		{
			if( it->second != rhs )
			{
				return ZErrorMessage(ZErrorFileFormatBroken, "Douplicated attribute : " + lhs, __FILE__, __LINE__, __FUNCTION__);
			}
		}
		else attribute_tray[lhs] = rhs;

		zerr = getdelim(iss2, lhs, "=");
	}
	return ZErrorMessage();
}


void I_ReadData::putParamToReadOrWrite(vector<RWParam_void> & tray) const
{
	tray.clear();
	tray.push_back( RWParam_void( this->putTagLabel() ) );
}


ZErrorMessageReadingFile I_ReadData::readFile(const string& filename, const string& file_label)
{
	ifstream ifs(filename.c_str());
    if(!ifs)
    {
       	return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileNotFound, __FILE__, __LINE__, __FUNCTION__));
    }

    return ZErrorMessageReadingFile(filename, this->readStream(ifs, file_label) );
}


ZErrorMessage I_ReadData::readContents(istream& is, const string& label0, string& s)
{
    string label;
    ZErrorMessage zerr;
    while ( true )
    {
    	zerr = readLabel(is, label);
    	if( zerr.putErrorType() != ZErrorNoError ) return zerr;
		if( label.empty() )
		{
			return ZErrorMessage(ZErrorFileFormatBroken, "The label is empty", __FILE__, __LINE__, __FUNCTION__);
		}

		if( label.at(0) != '!' && label.at(0) != '?' ) break;
    }
    if( label != label0 )
    {
		return ZErrorMessage(ZErrorFileFormatBroken, "The correct label is \"" + label0 + ", not \"" + label + "\"", __FILE__, __LINE__, __FUNCTION__);
    }

    return readToLabelEnd(is, label, s);
}


ZErrorMessage I_ReadData::readStream(istream& ifs, const string& file_label)
{
	string s;
	ZErrorMessage zerr = readContents(ifs, file_label, s);
    if( zerr.putErrorType() != ZErrorNoError )
    {
    	return zerr;
    }

    return this->readData(s);
}


ZErrorMessage I_ReadData::readData(const string& s)
{
    ZErrorMessage zerr;

	vector<RWParam_void> param_tray;
	this->putParamToReadOrWrite(param_tray);

    vector<Int4> set_tray(param_tray.size(), 0);
    istringstream iss(s);
    ostringstream oss;
    zerr = readData(iss, param_tray, set_tray, 0, 0, &oss);
ZLOG_INFO( oss.str() + "\n" );
	if( zerr.putErrorType() != ZErrorNoError ) return zerr;

	assert( param_tray.size() == set_tray.size() );
	
	vector<Int4>::const_iterator itset = set_tray.begin();
	for(vector<RWParam_void>::const_iterator it=param_tray.begin(); it<param_tray.end(); it++, itset++)
	{
		zerr = this->checkIfDataAreSet(*it, *itset);
    	if( zerr.putErrorType() != ZErrorNoError ) return zerr;
		zerr = this->checkData(*it);
    	if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	}

	return ZErrorMessage();
}



ZErrorMessage I_ReadData::readData(istream& iss0, vector<RWParam_void>& param_tray,
		vector<Int4>& set_tray,
		const Int4& outwidth, const Int4& stage, ostream* os)
{
    assert( param_tray.size() == set_tray.size() );

    const Int4 isize = param_tray.size();
	// Copy.
    multimap<string, pair<RWParam_void, const Int4> > param_map;
	for(Int4 i=stage; i<isize; i++)
	{
		param_map.insert( pair<string, pair<RWParam_void, Int4> >(param_tray[i].putLabel(),
							pair<RWParam_void, const Int4>(param_tray[i], i) ) );
	}
	assert( (Int4)param_map.size() + stage == isize );
	
	multimap<string, pair<RWParam_void, const Int4> >::iterator itmap;

	string label, s, str;
    map<string, string> attribute_tray;
    ZErrorMessage zerr, zerr_check = ZErrorMessage();
    while ( readLabelAll(iss0, label).putErrorType() == ZErrorNoError )
    {
		if( label.empty() )
		{
			return ZErrorMessage(ZErrorFileFormatBroken, "The label is empty", __FILE__, __LINE__, __FUNCTION__);
		}
		if( label.at(0) == '!' )
		{
	    	continue;
		}
		if( label.at(0) == '?' )
		{
	    	continue;
		}

		istringstream iss1(label);
    	iss1 >> label;

    	zerr = readAttribute(iss1, attribute_tray);
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;

    	zerr = readToLabelEnd(iss0, label, s);
        if( zerr.putErrorType() != ZErrorNoError )
        {
        	return zerr;
        }

    	pair< multimap<string, pair<RWParam_void, const Int4> >::iterator,
    			multimap<string, pair<RWParam_void, const Int4> >::iterator> it_pair;
    	it_pair = param_map.equal_range( label );
    	if( it_pair.second == it_pair.first )
    	{
    		continue;
    	}

    	for(itmap=it_pair.first; itmap!=it_pair.second; itmap++)
    	{
        	if( itmap->second.first.putProperty().putAttribute().first != "" )
        	{
            	const map<string, string>::const_iterator itstrmap = attribute_tray.find(itmap->second.first.putProperty().putAttribute().first);
            	if( itstrmap == attribute_tray.end()
            			|| itstrmap->second != itmap->second.first.putProperty().putAttribute().second )
            	{
            		if( itmap->second.first.putProperty().putAttribute().second != "" ) continue;
            		else
            		{
                		itmap->second.first.setAttribute( pair<string, string>(itstrmap->first, itstrmap->second) );
            		}
            	}
        	}
        	break;
    	}
    	if( itmap == it_pair.second ) continue;

		set_tray[itmap->second.second]++;
		const bool output_flag = ( os != NULL && set_tray[itmap->second.second] <= putOutputNumber(itmap->second.first.putProperty()) );
if( output_flag )
{
	*os << string(outwidth*2, ' ') + label;
	if( itmap->second.first.putProperty().putAttribute().first != "" )
	{
		const map<string, string>::const_iterator itstrmap = attribute_tray.find(itmap->second.first.putProperty().putAttribute().first);
		*os << " " << itstrmap->first << "=\"" << itstrmap->second << "\"";
	}
	*os << ": ";
}

		istringstream iss2(s);

		if( itmap->second.first.putType() == VOIDDATA )
		{
			const Int4 isize0 = param_tray.size();
			this->setData(itmap->second.first.putProperty(), param_tray);
if( output_flag )
{
*os << endl;
}

			assert( itmap->second.first.putProperty().putAttribute().first != "" || (Int4)param_tray.size() > isize0 );	// Check if parameters are inserted.
			set_tray.resize(set_tray.size()+(param_tray.size()-isize0), 0 );

			zerr = readData(iss2, param_tray, set_tray, outwidth+1, isize0, (output_flag?os:NULL));
		}
		else
		{
			zerr = removeLabel(iss2, str);
			if( zerr.putErrorType() == ZErrorNoError )
			{
				istringstream iss3(str);
				zerr = itmap->second.first.setData(iss3, (output_flag?os:NULL));
			}
if( output_flag )
{
*os << endl;
}
		}
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	}

    return ZErrorMessage();
}



ZErrorMessage I_ReadData::checkData(const RWParam_void& param) const
{
	if( param.putLabel() == this->putTagLabel() ) return ZErrorMessage();
	return param.checkData();
}

ZErrorMessage I_ReadData::checkIfDataAreSet(const RWParam_void& param,
		const Int4& set_tray) const
{
	const string Label = param.putProperty().putLabelWithAttribute();
	const string LongLabel = param.putLongLabel();
	const Int4 max_label_num = param.putMaxMultiNumber();
	const Int4 min_label_num = param.putMinMultiNumber();
	if( set_tray > max_label_num )
	{
    	return ZErrorMessage(ZErrorDouplicatedLabel, LongLabel + ": <"+Label+"> appears more than "+num2str<Int4>(max_label_num) + " times", __FILE__, __LINE__, __FUNCTION__);
	}

	if( set_tray < min_label_num )
	{
    	return ZErrorMessage(ZErrorLabelNotFound, LongLabel + ": <"+Label+"> appears less than "+num2str<Int4>(min_label_num) + " times", __FILE__, __LINE__, __FUNCTION__);
    }

	return ZErrorMessage();
}

ZErrorMessage I_ReadData::checkIfDataIsSet(const RWParam_void& param,
							 const Int4& set_tray,
							 const string& attribute_value,
							 const Int4& min_number,
							 const Int4& max_number)
{
	if( set_tray > max_number )
	{
		const string Label = param.putProperty().putLabel() + " " + param.putProperty().putAttribute().first + "=\"" + attribute_value + "\"";
		const string LongLabel = param.putLongLabel();
    	return ZErrorMessage(ZErrorDouplicatedLabel, LongLabel + ": <" + Label + "> appears more than "+num2str<Int4>(max_number) + " times", __FILE__, __LINE__, __FUNCTION__);
	}
	else if( set_tray < min_number )
	{
		const string Label = param.putProperty().putLabel() + " " + param.putProperty().putAttribute().first + "=\"" + attribute_value + "\"";
		const string LongLabel = param.putLongLabel();
    	return ZErrorMessage(ZErrorLabelNotFound, LongLabel + ": <"+Label+"> appears less than "+num2str<Int4>(min_number) + " times", __FILE__, __LINE__, __FUNCTION__);
	}

	return ZErrorMessage();
}


map<string, Int4> I_ReadData::MAP_TO_REPLACE_NUM_THREAD()
{
	map<string, Int4> replace_num_thread;
#ifdef _OPENMP
	const Int4 MAX_THREAD_NUM = omp_get_num_procs();
	replace_num_thread.insert( map<string, Int4>::value_type( "MAX", MAX_THREAD_NUM ) );
	replace_num_thread.insert( map<string, Int4>::value_type( "DEFAULT", max(1, MAX_THREAD_NUM - 1 ) ) );
#endif
	return replace_num_thread;
}

ZErrorMessage I_ReadData::REPLACE_NUM_THREAD(istream& is, Int4& arg)
{
	arg = 1;
#ifdef _OPENMP
	ZErrorMessage zerr = I_ReadData::setValue(MAP_TO_REPLACE_NUM_THREAD(), is, arg);
	if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	if( arg > omp_get_num_procs() ) arg = omp_get_num_procs();
#endif
	return ZErrorMessage();
};
