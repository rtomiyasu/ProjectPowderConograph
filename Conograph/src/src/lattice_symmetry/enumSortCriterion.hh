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
#ifndef ENUMSORTCRITERION_HH_
#define ENUMSORTCRITERION_HH_


enum eSortCriterion{ SCM = 0, SCWuM = 1, SCRevM = 2, SCSymM = 3, SCNN = 4 };

inline Int4 putNumberOfSortCriterion()
{
	static const Int4 SortCriterionNum = 5;
	return SortCriterionNum;
}

inline const string& putLabel(const eSortCriterion& i)
{
	static const size_t SortCriterionNum = 5;
	static const string label[SortCriterionNum] = {"M", "Mwu", "Mrev", "Msym", "NN" };
	return label[(size_t)i];
}


#endif /*ENUMSORTCRITERION_HH_*/
