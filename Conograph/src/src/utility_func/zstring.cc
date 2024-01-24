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
#include "../zerror_type/error_out.hh"
#include "zstring.hh"

// Check if the string is blank.
bool is_blank(const string& s)
{
	istringstream strstream(s);

	string str;
	strstream >> str;
	
	return strstream.fail();
}

// Split the argument str by the delimiter delim.
void split(const string& str, vector<string>& word, const char& delim)
{
	string::size_type index = str.find_first_of(delim);
	
	word.clear();
	word.push_back( str.substr(0, index) );
	
	string::size_type index0;
	while( index != string::npos ){
		index0 = index + 1;
		index = str.find_first_of(delim, index0);
		word.push_back( str.substr(index0, index-index0) );
	}
}

// Split the argument str by the delimiter delim.
bool split_term(istream& iss, vector<string>& word)
{
	word.clear();
	char c;
	if( !iss.get(c) ) return false;
	else word.push_back("");
	
	Int4 index = 0;
	if( c != '+' ) word.push_back( string(1,c) );
	
	while( iss.get(c) )
	{
		if( c == '+' )
		{
			if( word[index].empty() || word[index] == "-" ) return false;
			index++;
			word.resize(index+1);
		}
		else if( c == '-')
		{
			if( word[index].empty() || word[index] == "-" ) return false;
			index++;
			word.push_back("-");
		}
		else
		{
			word[index] += c;
		}
	}
	return true;
}

string getFileExtension(const string& inputFilePathString)
{
	int stringSize = inputFilePathString.size();
	const char* filePath = inputFilePathString.c_str();
	bool isFindDotCharacter = false;
	while(0 < stringSize){
		if(filePath[--stringSize] == '.'){
			isFindDotCharacter = true;
			break;
		}
	}
	string returnString;
	
	if(false != isFindDotCharacter){
		returnString = inputFilePathString.substr(stringSize+1, inputFilePathString.size()-stringSize);
	}
	
	return returnString;
}

void removeFileExtension(const string& fname0, string& fname)
{
    // Remove the file extension.
    fname = fname0;
    string::size_type index = fname.rfind(".");    
    if( index != string::npos ) fname = fname.substr(0, index);
}

// Read the sentence before delim and put it in ans.
ZErrorMessage getdelim(istream& ifs, string& ans, const string& delim)
{
	ans.clear();
	if( delim.empty() ) return ZErrorMessage(ZErrorArgument, "Delimiter is empty", __FILE__, __LINE__, __FUNCTION__);
	
	const int ilength = delim.length();
	int count = 0;
	char c;
	while( ifs.get(c) )
	{
		if( c == delim.at(count) )
		{
			count++;
			if( count >= ilength )
			{
				return ZErrorMessage();
			}
		}
		else
		{
			int i;
			for(i=count-1; i>=0; i--)
			{
				if( delim.at(i) == c )
				{
					if( i == 0 || delim.substr(0, i) == delim.substr(count-i, i) ) break;
				}
			}
			if( i >= 0 )
			{
				ans += delim.substr(0, count-i);
				count = i + 1;
			}
			else
			{
				ans += delim.substr(0, count) + c;
				count = 0;
			}
		}
	}
	return ZErrorMessage(ZErrorDelimiterNotFound, "The Delimiter \""+delim+"\" is not found", __FILE__, __LINE__, __FUNCTION__);
}
