#include "arg_helpers.hh"
#include <cstdlib>

/**
 * Returns an integer index which is unique for each supported type,
 * to be used in switch statements etc. Not depending on typeid allows
 * for determination of type_ID at compile time
 */
template<> const type_ID type<double>::ID = DOUBLE;
template<> const type_ID type<int>::ID    = INT;
template<> const type_ID type<char*>::ID  = STRING;
template<> const type_ID type<bool>::ID   = BOOL;

/**
 * String literal with the type name
 */
template<> const char* type<double>::name = "double";
template<> const char* type<int>::name    = "int";
template<> const char* type<char*>::name  = "string";
template<> const char* type<bool>::name   = "bool";




/**
 * Returns the preferred format character for the provided type 
 *
 * @return the type's preferred format character
 */
template <typename T> char type<T>::fmt() {

	// return appropriate type name
	switch (ID) {
		case DOUBLE: return 'g'; break;
		case INT:    return 'd'; break;
		case STRING:
		case BOOL:   return 's'; break;
	}

	// otherwise print an error and quit
	fprintf(stderr,"bad type template\n");
	exit(-1);
}

/**
 * Writes the format string used for this type's value,
 * given a specified field width. If specified width is 0,
 * no field width is placed in the format code
 *
 * @param[out] str string in which format code is written
 * @param[in] n value field width
 */
template <typename T> void type<T>::fmt(char *str,int n) {
	if(n==0) {
		sprintf(str,"%%%c",fmt());
	} else {
		sprintf(str,"%%%d%c",n,fmt());
	}
}

/**
 * A bool-specfic implementation for printing a type's value
 * to a C-string, which inserts the strings "true" or
 * "false" as appropriate
 *
 * @param[in] str C-string to print to
 * @param[in] n specified field width
 * @param[in] val the value to print
 */
template<> void type<bool>::sprint(char* fh,int n,bool val) {
	char format[10];
	fmt(format,n);
	sprintf(fh,format,val?"true":"false");
}

/**
 * Prints a type's value to a C-string
 *
 * @param[in] str C-string to print to
 * @param[in] n specified field width
 * @param[in] val the value to print
 */
template<typename T> void type<T>::sprint(char* fh,int n,T val) {
	char format[10];
	fmt(format,n);
	sprintf(fh,format,val);
}

/**
 * A bool-specfic implementation for printing a type's value
 * to a file stream, which inserts the strings "true" or
 * "false" as appropriate
 *
 * @param[in] fh files stream to print to
 * @param[in] n specified field width
 * @param[in] val the value to print
 */
template<> void type<bool>::fprint(FILE* fh,int n,bool val) {
	char format[10];
	fmt(format,n);
	fprintf(fh,format,val?"true":"false");
}

/**
 * Prints a type's value to a file stream
 *
 * @param[in] fh file stream to print to 
 * @param[in] n specified field width
 * @param[in] val the value to print
 */
template<typename T> void type<T>::fprint(FILE* fh,int n,T val) {
	char format[10];
	fmt(format,n);
	fprintf(fh,format,val);
}

/**
 * Parses a C-string as an double, returning the value in
 * the second argument. Returns true if successful.
 *
 * @param[in] str C-string to parse
 * @param[out] val the parsed value
 * @return whether the string was parsed correctly
 */
template<> bool type<double>::parse(char* str,double &val) {
	char *p;
	val = strtod(str,&p);
	return !(p==NULL||p[0]!='\0');
}

/**
 * Parses a C-string as an int, returning the value in
 * the second argument. Returns true if successful.
 *
 * @param[in] str C-string to parse
 * @param[out] val the parsed value
 * @return whether the string was parsed correctly
 */
template<> bool type<int>::parse(char* str,int &val) {
	char *p;
	val = static_cast<int>(strtol(str,&p,10));
	return !(p==NULL||p[0]!='\0');
}

/**
 * "Parses" a C-string as a string, returning the value in
 * the second argument. Trivially returns true.
 *
 * @param[in] str C-string to parse
 * @param[out] val the parsed value
 * @return whether the string was parsed correctly
 */
template<> bool type<char*>::parse(char* str,char* &val) {
	val=str;
	return true;
}

// explicit instantiation
template class type<double>;
template class type<int>;
template class type<bool>;
template class type<char*>;
