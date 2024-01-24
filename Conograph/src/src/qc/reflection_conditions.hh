#ifndef _REFLECTION_CONDITION_HH_
#define _REFLECTION_CONDITION_HH_

#include "../RietveldAnalysisTypes.hh"
#include <assert.h>

class BravaisType;

// Declaration of a struct to store information about reflection conditions
class DataReflectionConditions
{
public:
	enum eAxisOrder{ AXIS_ABC, AXIS_CAB, AXIS_BCA, AXIS_BAC, AXIS_CBA, AXIS_ACB };

private:
	string type;
	string str_conditions;
    bool (*is_not_extinct)(const Int4& h, const Int4& k, const Int4& l);
    eAxisOrder axis_order;

public:
	DataReflectionConditions(const string& arg1, const string& arg2, bool(*arg3)(const Int4& h, const Int4& k, const Int4& l), const eAxisOrder arg4 = AXIS_ABC)
    {
    	type = arg1;
    	str_conditions = arg2;
    	is_not_extinct = arg3;
    	axis_order = arg4;
    }

    string putShortStringType() const;
    inline const string& putStringType() const { return type; };
    inline const string& putStringConditions() const { return str_conditions; };
    inline bool isNotExtinct(const Int4& h, const Int4& k, const Int4& l) const
    {
    	if( axis_order == AXIS_ABC ) return (*is_not_extinct)(h, k, l);
    	if( axis_order == AXIS_CAB ) return (*is_not_extinct)(k, l, h);
    	if( axis_order == AXIS_BCA ) return (*is_not_extinct)(l, h, k);
    	if( axis_order == AXIS_BAC ) return (*is_not_extinct)(k, h, l);
    	if( axis_order == AXIS_CBA ) return (*is_not_extinct)(l, k, h);
    	if( axis_order == AXIS_ACB ) return (*is_not_extinct)(h, l, k);
    	assert(false);
    	return true;
    }
};

Int4 putNumberOfTypesOfSystematicAbsences(const BravaisType& type);

// Declaration of a function to get information about reflection conditions.
const DataReflectionConditions& putInformationOnReflectionConditions(const BravaisType& brav_type, const Int4& irc_type);

#endif
